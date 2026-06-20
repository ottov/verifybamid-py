"""Precompute a per-sample genotype matrix aligned to a marker panel (for CHIPMIX).

CHIPMIX needs each sample's own genotype at the panel markers. Calling
`extract_self_geno` re-reads the whole VCF once per CRAM -- wasteful across 100k runs
when many samples share one genotype VCF. This builds the genotype matrix ONCE: rows
aligned to panel.parquet marker order, one int8 column per sample (0 missing / 1 homref
/ 2 het / 3 homalt, matching verifyBamID setGenotype). Then per-CRAM CHIPMIX just reads
that sample's single column.

  build-chip --vcf cohort.vcf.gz --panel panel.parquet --out chip.parquet
  verifybamid ... --chip-matrix chip.parquet     # looks up the @RG SM column
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from cyvcf2 import VCF


def build(vcf_path, panel_path):
    m = pq.read_table(panel_path)
    chrom = m.column("chrom").to_pylist()
    pos = m.column("pos").to_pylist()
    n_markers = m.num_rows
    key_to_idx = {(c, p): i for i, (c, p) in enumerate(zip(chrom, pos))}

    vcf = VCF(vcf_path)
    samples = list(vcf.samples)
    # geno[marker, sample] int8: default 0 (missing) for markers absent from the VCF
    geno = np.zeros((n_markers, len(samples)), dtype=np.int8)

    t0 = time.time()
    seen = 0
    for v in vcf:
        idx = key_to_idx.get((v.CHROM, v.POS))
        if idx is None:
            continue
        arr = v.genotype.array()                  # (n_samples, >=2), missing = -1
        a = arr[:, 0]
        b = arr[:, 1]
        valid = a >= 0                            # first allele present
        bb = np.where(b < 0, 0, b)               # trailing-missing half-call -> REF
        alt = (a > 0).astype(np.int8) + (bb > 0).astype(np.int8)
        geno[idx] = np.where(valid, 1 + alt, 0)  # homref=1/het=2/homalt=3; missing=0
        seen += 1
        if seen % 100_000 == 0:
            print(f"  ...{seen:,} panel markers filled ({time.time()-t0:.0f}s)",
                  file=sys.stderr)

    cols = {s: pa.array(geno[:, j], pa.int8()) for j, s in enumerate(samples)}
    table = pa.table(cols)
    print(f"chip matrix: {n_markers:,} markers x {len(samples):,} samples "
          f"({seen:,} markers matched in VCF)", file=sys.stderr)
    return table


def load_sample(chip_matrix_path, sample_id):
    """Read one sample's per-marker genotype column from a chip matrix, or None."""
    schema = pq.read_schema(chip_matrix_path)
    if sample_id not in schema.names:
        print(f"[chip] sample {sample_id} not in {chip_matrix_path}; CHIPMIX=NA",
              file=sys.stderr)
        return None
    col = pq.read_table(chip_matrix_path, columns=[sample_id]).column(0)
    return np.asarray(col, dtype=np.int8)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Precompute per-sample genotype matrix for CHIPMIX.")
    ap.add_argument("--vcf", required=True, help="VCF with genotypes (cohort callset / array)")
    ap.add_argument("--panel", required=True, help="panel.parquet from build_panel")
    ap.add_argument("--out", required=True, help="output chip matrix parquet")
    a = ap.parse_args(argv)
    pq.write_table(build(a.vcf, a.panel), a.out)
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
