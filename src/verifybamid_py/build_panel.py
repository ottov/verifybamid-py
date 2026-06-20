"""Stage 0: build the marker panel from a population VCF.

Reproduces verifyBamID 1.1.3's marker-acceptance logic exactly (see
git/verifybamid/src/VerifyBamID.cpp:64-171 and VcfFile.cpp:846-1024), so the
resulting panel matches the ``#SNPS`` count of the original ``.selfSM`` output.

Marker acceptance (verifyBamID, in order):
  1. autosomal only (chr1..chr22)
  2. AF and callRate computed from genotypes (this VCF has no AC/AN/AF INFO)
  3. skip if AF < minAF        (one-sided; AF is the alt-allele frequency)
  4. skip if callRate < minCallRate
  5. skip if multiallelic      (len(ALT) > 1) -- checked last
  (no SNP-vs-indel filter; indels are kept as first-base markers)

Genotype counting (verifyBamID computeAlleleCounts + setSample semantics):
  * GT whose first allele is missing ("./.", ".", "./1", ...) -> whole genotype
    is dropped (contributes nothing to AC or AN).
  * a trailing-missing half-call ("1/.", "0/.") backfills the missing allele as
    REFERENCE, because verifyBamID parses it with atoi(".") == 0. So "1/." is a
    het and "0/." is hom-ref.
  * AC = number of alleles with index > 0; AN = number of counted alleles.
  * AF = AC / (AN + 1e-6); callRate = AN / 2 / n_samples.
"""

from __future__ import annotations

import argparse
import re
import sys
import time

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from cyvcf2 import VCF

_LEADING_INT = re.compile(r"^\s*[+-]?(\d+)")


def _atoi(s: str) -> int:
    """C atoi: leading optional sign + digits, 0 if none."""
    m = _LEADING_INT.match(s)
    return int(m.group(0)) if m else 0


def is_autosome(chrom: str) -> bool:
    """Replicate verifyBamID VcfHelper::isAutosome (VcfFile.cpp:75-87).

    Returns true for a bare positive int OR "chr"+positive-int -- which
    deliberately ACCEPTS chr23/24/25/26 (X/Y/XY/MT in numeric-coded VCFs),
    matching verifyBamID's behaviour. These contigs usually don't match the
    CRAM's chrX/Y/M names, so they carry no reads but still count in #SNPS.
    """
    if _atoi(chrom) > 0:
        return True
    if chrom[:3] == "chr" and _atoi(chrom[3:]) > 0:
        return True
    return False


def _genotype_array(variant):
    """Return an (n_samples, >=2) int array of allele indices, missing = -1.

    Uses cyvcf2's fast C path when available, falling back to the list API.
    """
    try:
        return variant.genotype.array()
    except AttributeError:
        return np.asarray(variant.genotypes, dtype=np.int32)


def allele_counts(gt: np.ndarray) -> tuple[int, int]:
    """(AC, AN) replicating verifyBamID computeAlleleCounts + setSample.

    gt[:, 0] / gt[:, 1] are the two allele indices; any negative value is a
    missing allele.
    """
    a = gt[:, 0]
    b = gt[:, 1]
    valid = a >= 0  # first allele present (GT does not start with '.')
    av = a[valid]
    bv = b[valid]
    bv = np.where(bv < 0, 0, bv)  # atoi('.') == 0 -> trailing missing -> REF
    an = int(2 * valid.sum())
    ac = int((av > 0).sum() + (bv > 0).sum())
    return ac, an


def build(
    vcf_path: str,
    min_af: float = 0.01,
    min_callrate: float = 0.50,
    snps_only: bool = False,
):
    vcf = VCF(vcf_path)
    n_samples = len(vcf.samples)
    if n_samples == 0:
        sys.exit("VCF has no genotypes; cannot derive allele frequencies.")

    chroms, poss, refs, alts, afs = [], [], [], [], []
    stats = dict(
        total=0, non_autosomal=0, multiallelic=0, no_af=0,
        af_fail=0, callrate_fail=0, kept=0,
    )
    t0 = time.time()
    for v in vcf:
        stats["total"] += 1
        if stats["total"] % 200_000 == 0:
            print(
                f"  ...{stats['total']:,} records  kept={stats['kept']:,}  "
                f"({time.time() - t0:.0f}s)",
                file=sys.stderr,
            )

        if not is_autosome(v.CHROM):
            stats["non_autosomal"] += 1
            continue

        ac, an = allele_counts(_genotype_array(v))
        if an == 0:
            stats["no_af"] += 1
            continue
        af = ac / (an + 1e-6)
        callrate = an / 2.0 / n_samples

        # verifyBamID applies AF then callRate BEFORE the biallelic check.
        if af < min_af:
            stats["af_fail"] += 1
            continue
        if callrate < min_callrate:
            stats["callrate_fail"] += 1
            continue
        if v.ALT is None or len(v.ALT) != 1:
            stats["multiallelic"] += 1
            continue
        if snps_only and not (len(v.REF) == 1 and len(v.ALT[0]) == 1):
            continue

        chroms.append(v.CHROM)
        poss.append(v.POS)                 # 1-based, matches VCF
        refs.append(v.REF[0])              # verifyBamID uses sRef[0]
        alts.append(v.ALT[0][0])           # ...and asAlts[0][0]
        afs.append(af)
        stats["kept"] += 1

    table = pa.table(
        {
            "chrom": pa.array(chroms, pa.string()),
            "pos": pa.array(poss, pa.int32()),
            "ref": pa.array(refs, pa.string()),
            "alt": pa.array(alts, pa.string()),
            "af": pa.array(afs, pa.float64()),
        }
    )
    return table, stats


def main(argv=None):
    p = argparse.ArgumentParser(description="Build verifyBamID marker panel from a VCF.")
    p.add_argument("--vcf", required=True)
    p.add_argument("--out", help="output parquet path (omit to only report counts)")
    p.add_argument("--min-af", type=float, default=0.01)
    p.add_argument("--min-callrate", type=float, default=0.50)
    p.add_argument("--snps-only", action="store_true",
                   help="restrict to biallelic SNPs (fast-mode panel; not verifyBamID-exact)")
    p.add_argument("--expect", type=int, default=None,
                   help="assert kept marker count equals this (validation gate)")
    args = p.parse_args(argv)

    table, stats = build(args.vcf, args.min_af, args.min_callrate, args.snps_only)

    print("\nmarker breakdown:", file=sys.stderr)
    for k, val in stats.items():
        print(f"  {k:14s} {val:,}", file=sys.stderr)

    if args.out:
        pq.write_table(table, args.out)
        print(f"\nwrote {stats['kept']:,} markers -> {args.out}", file=sys.stderr)

    if args.expect is not None:
        if stats["kept"] != args.expect:
            sys.exit(f"FAIL: kept {stats['kept']:,} != expected {args.expect:,}")
        print(f"PASS: kept == {args.expect:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
