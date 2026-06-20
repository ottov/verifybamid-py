"""End-to-end driver: CRAM (local / s3:// / presigned https) -> .selfSM.

Streams the CRAM (Stage 1 pileup) and estimates FREEMIX/FREELK1/FREELK0 (Stage 2),
writing a verifyBamID-compatible .selfSM row. For s3:// inputs it presigns the CRAM
and its .crai so the pysam-bundled libcurl can range-read them directly.

Example:
  verifybamid --cram s3://bucket/sample.cram --ref hg38.fa \
              --panel panel.parquet --out results/sample --jobs 18
"""

from __future__ import annotations

import argparse
import os
import sys
from urllib.parse import urlparse

from . import estimate, pileup

SELF_SM_HEADER = [
    "#SEQ_ID", "RG", "CHIP_ID", "#SNPS", "#READS", "AVG_DP", "FREEMIX", "FREELK1",
    "FREELK0", "FREE_RH", "FREE_RA", "CHIPMIX", "CHIPLK1", "CHIPLK0", "CHIP_RH",
    "CHIP_RA", "DPREF", "RDPHET", "RDPALT",
]


def _presign(s3_uri: str, expires: int, region: str) -> str:
    import boto3
    u = urlparse(s3_uri)
    client = boto3.client("s3", region_name=region)
    return client.generate_presigned_url(
        "get_object", Params={"Bucket": u.netloc, "Key": u.path.lstrip("/")},
        ExpiresIn=expires)


def resolve_urls(cram: str, crai: str | None, expires: int, region: str):
    """Return (cram_url, crai_url) ready for pysam, presigning s3:// inputs."""
    if cram.startswith("s3://"):
        crai_src = crai or (cram + ".crai")
        return _presign(cram, expires, region), _presign(crai_src, expires, region)
    return cram, (crai or (cram + ".crai"))


def sample_id_from_cram(cram_url: str, ref: str, crai_url: str | None) -> str:
    af = pileup.open_cram(cram_url, ref, crai_url)
    try:
        for rg in af.header.get("RG", []):
            if rg.get("SM"):
                return rg["SM"]
    finally:
        af.close()
    base = os.path.basename(urlparse(cram_url).path or cram_url)
    return base.split(".cram")[0]


def main(argv=None):
    p = argparse.ArgumentParser(description="CRAM -> verifyBamID .selfSM (FREEMIX/FREELK).")
    p.add_argument("--cram", required=True, help="local path, s3:// URI, or presigned https URL")
    p.add_argument("--crai", help="explicit .crai (default: <cram>.crai)")
    p.add_argument("--ref", required=True)
    p.add_argument("--panel", required=True)
    p.add_argument("--out", required=True, help="output prefix (writes <out>.selfSM)")
    p.add_argument("--seq-id", help="sample id (default: @RG SM, else filename)")
    p.add_argument("--chip-vcf", help="VCF with the sample's genotypes -> also compute CHIPMIX")
    p.add_argument("--chip-matrix", help="precomputed chip matrix (build-chip) -> CHIPMIX")
    p.add_argument("--chip-id", help="sample id within chip source (default: --seq-id)")
    p.add_argument("--best", action="store_true",
                   help="(optional) identity/sample-swap check: best-matching individual "
                        "in the chip matrix; requires --chip-matrix")
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--contig", action="append", dest="contigs",
                   help="restrict to contig(s); repeatable (mainly for testing)")
    p.add_argument("--max-span", type=int,
                   help="cap fetch span (bp); set ~150000 with a downsampled panel")
    p.add_argument("--keep-pileup", action="store_true",
                   help="keep <out>.pileup/.markers.parquet (default: write then reuse)")
    p.add_argument("--region", default=os.environ.get("AWS_REGION", "us-east-1"))
    p.add_argument("--expires", type=int, default=3600)
    p.add_argument("--grid", type=float, default=0.05)
    p.add_argument("--min-mapq", type=int)
    p.add_argument("--min-q", type=int)
    p.add_argument("--max-q", type=int, default=40)
    p.add_argument("--max-depth", type=int)
    a = p.parse_args(argv)

    cram_url, crai_url = resolve_urls(a.cram, a.crai, a.expires, a.region)
    seq_id = a.seq_id or sample_id_from_cram(cram_url, a.ref, crai_url)
    print(f"[verifybamid] sample={seq_id}  source={a.cram}", file=sys.stderr)

    # Stage 1: stream + pile
    counts = pileup.run(cram_url, a.ref, a.panel, crai_url=crai_url, out=a.out,
                        jobs=a.jobs, contigs=a.contigs, max_span=a.max_span,
                        min_mapq=a.min_mapq, min_q=a.min_q, max_q=a.max_q,
                        max_depth=a.max_depth)

    # Stage 2: estimate
    d = estimate.load(f"{a.out}.markers.parquet", f"{a.out}.pileup.parquet", max_q=a.max_q)
    freemix, freelk1, freelk0 = estimate.optimize(d, d["gfo"], grid=a.grid, max_alpha=0.5)

    chip = ["NA", "NA", "NA"]
    chip_id_out = "NA"
    if a.best:
        if not a.chip_matrix:
            sys.exit("--best requires --chip-matrix")
        best_id, cmix, clk1, clk0 = estimate.optimize_best(d, a.chip_matrix, grid=a.grid)
        chip = [f"{cmix:.5f}", f"{clk1:.2f}", f"{clk0:.2f}"]
        chip_id_out = best_id
        if best_id != seq_id:
            print(f"[verifybamid] ** SWAP: best-match {best_id} != seq-id {seq_id} **",
                  file=sys.stderr)
    elif a.chip_matrix or a.chip_vcf:
        if a.chip_matrix:
            from . import build_chip
            sg = build_chip.load_sample(a.chip_matrix, a.chip_id or seq_id)
        else:
            sg = estimate.extract_self_geno(a.chip_vcf, a.chip_id or seq_id,
                                            f"{a.out}.markers.parquet")
        if sg is not None:                       # None = sample absent -> CHIPMIX NA
            cmix, clk1, clk0 = estimate.optimize(d, estimate.chip_weights(d, sg),
                                                 grid=a.grid, max_alpha=0.95)
            chip = [f"{cmix:.5f}", f"{clk1:.2f}", f"{clk0:.2f}"]
            chip_id_out = a.chip_id or seq_id

    row = [seq_id, "ALL", chip_id_out, str(counts["n_snps"]), str(counts["n_reads"]),
           f"{counts['avg_dp']:.2f}", f"{freemix:.5f}", f"{freelk1:.2f}",
           f"{freelk0:.2f}", "NA", "NA", chip[0], chip[1], chip[2],
           "NA", "NA", "NA", "NA", "NA"]
    self_sm = f"{a.out}.selfSM"
    with open(self_sm, "w") as fh:
        fh.write("\t".join(SELF_SM_HEADER) + "\n")
        fh.write("\t".join(row) + "\n")

    if not a.keep_pileup:
        for suffix in (".pileup.parquet", ".markers.parquet"):
            try:
                os.remove(f"{a.out}{suffix}")
            except OSError:
                pass

    chip_msg = f"  CHIPMIX={chip[0]}" if a.chip_vcf else ""
    print(f"[verifybamid] FREEMIX={freemix:.5f}  FREELK1={freelk1:.2f}  "
          f"FREELK0={freelk0:.2f}{chip_msg}  -> {self_sm}", file=sys.stderr)
    print("\t".join(SELF_SM_HEADER))
    print("\t".join(row))


if __name__ == "__main__":
    main()
