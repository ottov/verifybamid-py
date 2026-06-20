"""Build a sparse 'fast-mode' panel from a full panel.parquet.

FREEMIX is a single contamination fraction and is statistically saturated by a few
thousand common, well-spaced SNPs, so we can drop from ~484k markers to ~10-20k with
negligible effect on the estimate. The point is I/O: with markers spaced wider than a
CRAM slice (~50-80 kb), targeted .crai fetches skip most slices instead of reading the
whole CRAM. (FREELK values won't match the full-panel run -- different marker set.)

Selection:
  * autosomes chr1-22 only (chr23-26 carry no reads anyway)
  * common SNPs: MAF = min(af, 1-af) >= --min-maf  (informative for contamination)
  * then N markers taken evenly across the position-sorted set (genome-wide spread)
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

AUTOSOMES = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}


def downsample(panel_path, n, min_maf):
    t = pq.read_table(panel_path)
    chrom = np.array(t.column("chrom").to_pylist())
    af = np.asarray(t.column("af"), dtype=np.float64)
    maf = np.minimum(af, 1 - af)

    keep = np.array([c in AUTOSOMES for c in chrom]) & (maf >= min_maf)
    idx = np.flatnonzero(keep)  # already position-sorted within contig
    if len(idx) > n:
        # even spacing across the sorted common-SNP set -> ~uniform genomic spread
        pick = np.unique(np.linspace(0, len(idx) - 1, n).round().astype(int))
        idx = idx[pick]
    return t.take(pa.array(idx)), len(idx)


def main(argv=None):
    p = argparse.ArgumentParser(description="Downsample a full panel to a sparse fast-mode panel.")
    p.add_argument("--panel", required=True, help="full panel.parquet from build_panel")
    p.add_argument("--out", required=True, help="output sparse panel.parquet")
    p.add_argument("-n", "--num", type=int, default=20000, help="target marker count")
    p.add_argument("--min-maf", type=float, default=0.10)
    a = p.parse_args(argv)

    tbl, kept = downsample(a.panel, a.num, a.min_maf)
    pq.write_table(tbl, a.out)
    pos = np.asarray(tbl.column("pos"))
    span = int(pos.max() - pos.min()) if len(pos) else 0
    print(f"kept {kept:,} markers (MAF>={a.min_maf}, autosomal) -> {a.out}", file=sys.stderr)
    print(f"~1 marker per {3_100_000_000 // max(kept,1):,} bp (slices ~50-80kb)", file=sys.stderr)


if __name__ == "__main__":
    main()
