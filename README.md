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
submits one streaming SLURM job per sample. Needs **Java 17+,
[Nextflow](https://www.nextflow.io/docs/latest/install.html), [`uv`](https://docs.astral.sh/uv/),
and the `aws` CLI**, plus a **shared filesystem** that every compute node can see.

### 1. Install (once, on the shared filesystem)

```bash
git clone https://github.com/ottov/verifybamid-py.git
cd verifybamid-py
uv sync                       # builds .venv/ with the verifybamid commands
```

If other users will run from your install, keep uv's Python inside the checkout.
Otherwise `.venv/bin/python` points into your `~/.local/share/uv`, and their jobs fail
with `bad interpreter: Permission denied` (exit 126):

```bash
UV_PYTHON_INSTALL_DIR=$PWD/.uv-python uv sync
chmod -R a+rX .
```

### 2. Check the prerequisites

- The reference FASTA the CRAMs were aligned to, with its `.fai`, on the shared filesystem.
- AWS credentials that can read the CRAMs and the panel VCF **from the compute nodes**
  (e.g. an NFS-shared `~/.aws`). Check with `aws s3 ls s3://bucket/path/sample.cram`.

### 3. Write the sample sheet

Tab-separated, no header, one sample per line:

```
sample_id <TAB> s3://bucket/project/PANEL_PREFIX <TAB> s3://bucket/path/sample.cram
```

`PANEL_PREFIX` is the project's population VCF **without** `.vcf.gz`; the pipeline
appends it. For `s3://bucket/project/cohort.vcf.gz`, column 2 is
`s3://bucket/project/cohort`. Every sample with the same prefix shares one panel.

### 4. Point it at your cluster

The bundled config submits to a partition named `defaultq`. Put your cluster's settings
in a small `site.config`; any config file given later on the command line wins:

```groovy
process {
    queue          = 'your_partition'
    clusterOptions = '--qos=your_qos'                    // any extra sbatch flags
    beforeScript   = 'export PATH=/usr/local/bin:$PATH'  // if aws isn't on the jobs' PATH
}
```

### 5. Launch

Launch from a project directory on the shared filesystem. `work/`, `panels/` and
`results/` are created in the directory you launch from. Use absolute paths to the
checkout. The Nextflow driver has to run for the whole batch, so start it with `nohup`
on a login node, or submit it as its own long, 1-CPU job. Don't start it inside an
interactive session that will end.

```bash
VBID=/shared/path/to/verifybamid-py
cd /shared/path/to/project

nohup nextflow run $VBID/nextflow/main.nf \
  -c $VBID/nextflow/nextflow.config -c site.config \
  --samples samples.tsv \
  --ref     /shared/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --bindir  $VBID/.venv/bin \
  -resume > run.log 2>&1 &
```

`--ref` and `--bindir` are required: the defaults in `nextflow.config` are one
developer's paths.

The run builds the marker panel and CHIP matrix **once per panel** (cached in `panels/`
and reused on re-runs), then streams each CRAM as an independent SLURM job. If the
driver stops for any reason, run the same command again; `-resume` skips finished samples.

Outputs land in `results/`:

- `contamination.selfSM` — merged table, one row per sample (`FREEMIX`, `FREELK1/0`, `CHIPMIX`, …)
- `<sample>.selfSM` — one file per sample
- `failed_samples.txt` — samples that failed after retries (empty if all passed; re-run by
  trimming the sheet to these and re-launching with `-resume`)
- `_report.html` / `_timeline.html` / `_trace.txt` — run stats (`_trace.txt` has per-sample durations)

### Knobs

Pass as `--flag value` on the command line, or set them in `site.config` under `params { }`.
`max_streams` is read at launch, so changing it means stopping the driver and relaunching
with `-resume`. Samples in progress when you stop it are cancelled and re-run.

| flag | default | when to change |
|---|---|---|
| `--max_streams` | `8` | CRAMs streaming at once, which is how many samples run at a time. Raise it to go faster; see "Throughput and egress". |
| `--cpus` | `2` | cores (parallel fetches) per sample; 2–4 is the sweet spot. |
| `--fast` | `true` | `false` uses the full panel (the same markers as stock verifyBamID): slower, more S3 reads. |
| `--fast_n` | `20000` | markers in the fast panel; 20k is well-validated. |
| `--chipmix` | `true` | `false` skips the CHIP matrix. Samples missing from the panel VCF get `CHIPMIX=NA` either way. |
| `--region` | `us-east-1` | the S3 buckets' region. |
| `-profile local` | (slurm) | run on one node instead of submitting to SLURM (testing). |

### Troubleshooting

| symptom | cause |
|---|---|
| `Remote resource not found: https://api.github.com/repos/.../main.nf` | The path to `main.nf` doesn't exist from where you launched, so Nextflow looked for a GitHub project with that name. Use the absolute path. |
| `sbatch` rejects the partition | The bundled `defaultq` doesn't exist on your cluster. Set `process.queue` in `site.config`. |
| `aws: command not found` in `PANEL` / `CHIP` | Set `beforeScript` in `site.config`. |
| exit 126, `bad interpreter: Permission denied` | Another user's jobs can't read your uv Python. See step 1. |
| `Cannot compare java.lang.String with value '16'` | Your checkout predates the fix for Nextflow 26+, where command-line values arrive as text. On those checkouts `--fast false` is also silently ignored. Run `git pull`, or set the values in `site.config` under `params { }`. Mid-batch, use the config: the pull changes the per-sample step, so `-resume` would re-run samples that already finished. |
| sample fails with `likely a degraded S3 stream` | Its download under-read (coverage guard). Retried automatically; if it's frequent, lower `--max_streams`. |

## Install (single-sample use)

```bash
uv sync
```

Requires a local reference FASTA (same one the CRAMs were aligned to) to decode CRAM.
The rest of this README covers the single-sample / per-stage commands the pipeline calls.

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
apart genome-wide and still touch most CRAM slices (see "Throughput and egress" below):

```bash
uv run downsample --panel panel.parquet --out fast20k.parquet -n 20000 --min-maf 0.10
```

### 2+3. Estimate (end-to-end)

```bash
uv run verifybamid \
  --cram  s3://bucket/sample.cram \
  --ref   GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --panel panel.parquet \
  --out   results/sample \
  --jobs  4
```

- `--cram` takes a local path, an `s3://` URI, or a presigned https URL.
- Writes `results/sample.selfSM`.
- `--jobs` is parallel fetches; streaming is bandwidth-bound, not CPU-bound.
- With a downsampled panel, add `--max-span 1000000` (the measured egress optimum).
- Add `--chip-vcf cohort.vcf.gz` to also compute CHIPMIX if the sample is in it.

For `s3://` inputs the CRAM and its `.crai` are presigned via boto3 (region from
`--region` / `AWS_REGION`, default `us-east-1`). Run **in-region** for free, fast egress.

### Estimate from a VCPA marker table

The VCPA pipeline's marker pileup already records, per panel marker, the depth and
phred sum of each allele. `estimate --table` reads that directly -- no CRAM access,
~5 s per sample:

```bash
uv run estimate --markers markers.parquet --seq-id SM \
  --table SM.ad.tsv.gz                     # full table
  # or  --table SM.markers.manifest.tsv    # shard set, via its manifest
  # or  --table path/to/markers/           # directory holding one manifest
  # or  --table s3://bucket/.../markers/   # shard set on S3 (prefix or manifest key)
```

Table rows are `chrom pos ref alt refDepth altDepth refQsum altQsum`. Each allele is
scored at its mean quality `round(Qsum/depth)`, weighted by its depth, which gives
the same likelihood as expanding it to one pileup row per read. A shard set is read
only through its manifest and must match it exactly: every shard present with the
recorded bytes and rows, and the manifest covering the whole panel. Otherwise
`estimate` exits non-zero without a result. So does a table covering fewer than
`--min-cov-frac` (default 0.80) of panel markers, since an empty pileup still
produces a plausible-looking FREEMIX. `#READS` and `AVG_DP` are filled from the table.

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
| `build-chip` | genotype VCF + panel → per-sample genotype matrix for CHIPMIX |
| `pileup` | stream CRAM → per-marker base/quality pileup |
| `estimate` | pileup or VCPA marker table → FREEMIX/FREELK (+CHIPMIX) |
| `verifybamid` | end-to-end: CRAM → .selfSM |

## Throughput and egress

Measured streaming from S3 to on-premises SLURM nodes (fast 20k panel, `max_span=1M`,
bytes counted off the NIC):

| metric | value | note |
|---|---|---|
| egress / sample | **~14.3 GB** | of a ~19 GB CRAM → ~25% saving, not order-of-magnitude |
| `max_span` optimum | **1,000,000** | tighter re-downloads shared CRAM slices; uncapped pulls marker-free gaps |
| `fast_n` lever | sub-linear | 20k→5k markers saves only ~30% egress, for real accuracy loss |
| one stream | **~70 Mbit/s** | ~20–30 min per sample at `--cpus 2` |

Throughput is streams × per-stream rate. It stops scaling when one of two links fills:

- **Each node's network link.** Slurm packs samples onto nodes: at `--cpus 2`, an 8-CPU
  node takes 4 of them (~280 Mbit/s). A 1 Gbit/s node tops out around 600 Mbit/s in
  practice, so 10 streams on one node each ran several times slower. Check placement
  with `squeue -u $USER -o %N`.
- **The site's path to S3.** This is shared by every node. On-premises, it was not the
  limit at ~900 Mbit/s across four nodes: the streams already running kept their full
  rate. Its real ceiling depends on your site.

To tune `--max_streams`, raise it in steps (8 → 16 → 24), relaunching with `-resume`
each time. After each step, check the per-sample durations in `results/_trace.txt`. If
they stay flat, the extra streams are adding throughput. If they rise, or samples start
failing the coverage guard, you've hit a ceiling; go back a step.

Egress is the same whatever the concurrency; only the finish date changes. Projected for
**100k samples**:

- egress ≈ 14.3 GB × 100k ≈ **1.4 PB**, billed as S3 internet egress when pulled out of region
- transfer time ≈ 1.4 PB ÷ aggregate bandwidth: **~7 months at 600 Mbit/s, ~4.4 months
  at 1 Gbit/s, ~2 weeks at 10 Gbit/s**

Practical levers, biggest first:

1. **Run compute in-region (AWS).** Egress becomes free and each node gets multi-Gbps to
   S3, so months become days. This is by far the largest lever.
2. **More aggregate bandwidth on-prem:** spread streams across more nodes (a few per
   1 Gbit/s node) until the site's path to S3 fills; beyond that, a fatter S3↔site link
   (Direct Connect / more WAN) helps roughly linearly.
3. Fewer markers / a spatially-clustered panel would cut egress further, but the first
   trades accuracy and the second needs validation that the contamination model holds.

## Notes

- Deliverable columns: `FREEMIX`, `FREELK1`, `FREELK0`, and `CHIPMIX`/`CHIPLK*` when
  genotypes are available; reference-bias columns are fixed (`--free-mix` mode) and emit `NA`.
- Behavior follows the verifyBamID 1.1.3 C++ source
  ([statgen/verifyBamID](https://github.com/statgen/verifyBamID)), which was read to
  match it exactly.
