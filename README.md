# verifybamid-py

Streaming reimplementation of [verifyBamID](https://genome.sph.umich.edu/wiki/VerifyBamID)
1.1.3's chip-free contamination estimation (`FREEMIX` / `FREELK1` / `FREELK0`, plus
`CHIPMIX` when sample genotypes are available), built to **stream CRAMs directly from
S3** instead of staging the whole file to disk — for contamination QC across ~100k WGS
CRAMs.

Validated bit-for-bit against verifyBamID 1.1.3 across the full contamination spectrum:
clean (FREEMIX 0.000), low (0.014), and high (0.454), with FREELK matching to the cent
on the matching panel and CHIPMIX to 1e-5.

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
fetches skip slices — ~20x less data with the contamination call preserved:

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
  --jobs 18 \
  [--max-span 150000] \                 # set with a downsampled panel (targeted fetch)
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

## Notes

- Deliverable columns: `FREEMIX`, `FREELK1`, `FREELK0`, and `CHIPMIX`/`CHIPLK*` when
  genotypes are available; reference-bias columns are fixed (`--free-mix` mode) and emit `NA`.
- The reference C++ source (read to match behavior exactly) lives at `../git/verifybamid`.
