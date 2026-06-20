# verifybamid-py

Streaming reimplementation of [verifyBamID](https://genome.sph.umich.edu/wiki/VerifyBamID)
1.1.3's chip-free contamination estimation (`FREEMIX` / `FREELK1` / `FREELK0`, plus
`CHIPMIX` when sample genotypes are available), built to **stream CRAMs directly from
S3** instead of staging the whole file to disk — for contamination QC across ~100k WGS
CRAMs.

Validated bit-for-bit against verifyBamID 1.1.3 across the full contamination spectrum:
clean (FREEMIX 0.000), low (0.014), and high (0.454), with FREELK matching to the cent
on the matching panel and CHIPMIX to 1e-5.

## Cluster quickstart (SLURM)

Run a whole project (10,000s of CRAMs sharing one panel) as a Nextflow pipeline that
submits one streaming SLURM job per sample. Needs **Java 17+, Nextflow, `uv`, and the
`aws` CLI** on the cluster, plus a **shared filesystem** visible to all compute nodes.

```bash
# 1. clone + install
git clone git@github.com:ottov/verifybamid-py.git
cd verifybamid-py
uv sync                       # builds .venv/ with the verifybamid commands

# 2. one-time prerequisites
ls /shared/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa{,.fai}   # ref + .fai, shared
aws configure                 # AWS creds readable on COMPUTE nodes (NFS-shared ~/.aws)

# 3. sample sheet: TSV, NO header, tab-separated
#    sample_id <TAB> s3://bucket/project/PANEL_PREFIX <TAB> s3://bucket/path/sample.cram
#    (".vcf.gz" is appended to PANEL_PREFIX; the panel is shared by all samples in it)

# 4. launch (submits to SLURM by default); run from the repo root
nextflow run nextflow/main.nf -c nextflow/nextflow.config \
  --samples /path/to/sample_list.tsv \
  --ref     /shared/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --bindir  "$PWD/.venv/bin" \
  --region  us-east-1 \
  -resume
```

The run builds the marker panel + CHIP matrix **once per project** (cached in `panels/`,
reused on re-runs) and streams each CRAM as an independent SLURM job. Outputs land in
`results/`:

- `contamination.selfSM` — merged table, one row per sample (`FREEMIX`, `FREELK1/0`, `CHIPMIX`, …)
- `failed_samples.txt` — samples that failed after retries (empty if all passed; re-run by
  trimming the sheet to these and re-launching with `-resume`)
- `_report.html` / `_timeline.html` / `_trace.txt` — run stats

Knobs (as `--flag value`, or edit `nextflow/nextflow.config`):

| flag | default | when to change |
|---|---|---|
| `--max_streams` | `8` | concurrent S3 streams. **Raise a lot in-region/AWS**; lower if coverage-guard retries appear over a thin WAN. |
| `--cpus` | `2` | cores per sample; 2–4 is the sweet spot. |
| `--chipmix` | `true` | set `false` if the sample isn't in the project genotype VCF (then `FREEMIX` only). |
| `--fast_n` | `20000` | markers; 20k is well-validated. |
| `-profile local` | (slurm) | run on one node instead of submitting to SLURM (testing). |

> **It's WAN-bandwidth-bound** (~14 GB/sample), not CPU-bound — see "Scaling / egress"
> below. `-resume` makes failures/preemptions cheap to recover.

The rest of this README covers the single-sample / per-stage commands the pipeline calls.

## Install

```bash
uv sync
```

Requires a local reference FASTA (same one the CRAMs were aligned to) to decode CRAM.

## Pipeline

Three stages, exposed as one end-to-end command plus per-stage commands.

### 1. Build the marker panel (once per population VCF)

```bash
uv run build-panel --vcf panel.vcf.gz --out panel.parquet
```

Replicates verifyBamID's marker acceptance (autosomal incl. chr23-26 numeric coding,
biallelic, AF>=0.01 one-sided, callRate>=0.50, AF from genotypes via the exact
`computeAlleleCounts`/`setSample` semantics). The panel is shared across all samples.

Optional sparse "fast-mode" panel (≈20k common, well-spaced SNPs) so targeted `.crai`
fetches skip slices, with the contamination call preserved. The saving is real but
modest — **measured ~14 GB of a ~19 GB CRAM (~25%)**, because 20k markers sit ~155 kb
apart genome-wide and still touch most CRAM slices (see "Scaling / egress" below):

```bash
uv run downsample --panel panel.parquet --out fast20k.parquet -n 20000 --min-maf 0.10
```

### 2+3. Estimate (end-to-end)

```bash
uv run verifybamid \
  --cram s3://bucket/sample.cram \      # local path, s3:// URI, or presigned https URL
  --ref  GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --panel panel.parquet \
  --out  results/sample \               # writes results/sample.selfSM
  --jobs 4 \                            # streaming is WAN-bandwidth-bound, not CPU-bound
  [--max-span 1000000] \                # with a downsampled panel; 1M is the egress optimum
  [--chip-vcf cohort.vcf.gz]            # also compute CHIPMIX if the sample is in it
```

For `s3://` inputs the CRAM and its `.crai` are presigned via boto3 (region from
`--region` / `AWS_REGION`, default `us-east-1`). Run **in-region** for free, fast egress.

## CHIPMIX

CHIPMIX needs the sample's own genotypes (it compares reads against a known genotype).
Pass `--chip-vcf` pointing at any VCF containing the sample (its cohort callset or an
external SNP array; decoupled from `--panel`). If the sample isn't there, `CHIPMIX=NA`
and `FREEMIX` (which needs no per-sample genotypes) is the contamination estimate.

## Commands

| command | purpose |
|---|---|
| `build-panel` | population VCF → marker panel (chrom,pos,ref,alt,af) |
| `downsample` | full panel → sparse fast-mode panel |
| `pileup` | stream CRAM → per-marker base/quality pileup |
| `estimate` | pileup → FREEMIX/FREELK (+CHIPMIX) |
| `verifybamid` | end-to-end: CRAM → .selfSM |

## Scaling / egress

Measured on a real HPC node streaming from S3 over the WAN (fast 20k panel,
`max_span=1M`, bytes counted off the NIC):

| metric | value | note |
|---|---|---|
| egress / sample | **~14.3 GB** | of a ~19 GB CRAM → ~25% saving, not order-of-magnitude |
| `max_span` optimum | **1,000,000** | tighter re-downloads shared CRAM slices; uncapped pulls marker-free gaps |
| `fast_n` lever | sub-linear | 20k→5k markers saves only ~30% egress, for real accuracy loss |
| WAN aggregate | **~600 Mbit/s** (~75 MB/s) | the hard ceiling; saturated by a few concurrent streams |

**The binding constraint at scale is WAN bandwidth, not CPU.** Per-sample wall time and
core count barely matter: 10 concurrent streams on one node each slow ~7× because they
share the same ~75 MB/s pipe. The coverage guard catches streams that silently
under-read when the link is oversubscribed.

Projected for **100k samples** streaming on-prem:

- egress ≈ 14.3 GB × 100k ≈ **1.4 PB**
- transfer time ≈ 1.4 PB ÷ 75 MB/s ≈ **~7 months of continuous WAN transfer**, *regardless
  of how many cores or nodes you throw at it* — it's bandwidth-bound
- plus S3 internet-egress cost on ~1.4 PB if not pulled in-region

Practical levers, biggest first:

1. **Run compute in-region (AWS).** Egress becomes free and the pipe becomes multi-Gbps;
   the 7-month WAN bottleneck collapses to days. This is by far the largest lever.
2. **Fatter S3↔HPC link** (Direct Connect / more WAN) — improvement is ~linear in Gbps.
3. **Keep per-node concurrency low** (2–4 streams) so you stay under the throttle/coverage
   cliff; adding more just inflates per-sample latency without raising throughput.
4. Fewer markers / a spatially-clustered panel would cut egress further, but the first
   trades accuracy and the second needs validation that the contamination model holds.

## Notes

- Deliverable columns: `FREEMIX`, `FREELK1`, `FREELK0`, and `CHIPMIX`/`CHIPLK*` when
  genotypes are available; reference-bias columns are fixed (`--free-mix` mode) and emit `NA`.
- The reference C++ source (read to match behavior exactly) lives at `../git/verifybamid`.
