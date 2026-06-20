"""Stage 2: estimate FREEMIX/FREELK (and optionally CHIPMIX/CHIPLK) from the pileup.

Replicates verifyBamID 1.1.3:
  * FREE  (--free-mix, default): computeMixLLKs (VerifyBamID.cpp:510-641). The intended
    sample's genotype is marginalized over the HWE prior, same as the contaminant.
  * CHIP  (--chip-mix / --self):  computeIBDLLKs (VerifyBamID.cpp:646-803). The intended
    sample's genotype is FIXED to its known (chip/self-VCF) genotype; only the
    contaminant is marginalized over HWE.

Both share the per-marker base likelihood; they differ only in the final
marginalization weight. Per marker i with af clamped to [0.001,0.999] and observed
bases b_j (error e_j = 10^(-min(q_j,maxQ)/10)):

  pSN[ref]=[1,0.5,0]  pSN[alt]=[0,0.5,1]                 (P(allele|genotype 0/1/2))
  A[k1,k2]=fMix*pSN_ref[k1]+(1-fMix)*pSN_ref[k2]         (P draw ref allele)
  B[k1,k2]=fMix*pSN_alt[k1]+(1-fMix)*pSN_alt[k2]
  baseLK_j[k1,k2]=A*P(b_j|ref)+B*P(b_j|alt),  P(b|x)=(1-e) if b==x else e/3
  seg[k1,k2]=sum_j log baseLK_j                          (markerLK in log space)
  perMarker = sum_{k1,k2} exp(seg)[k1,k2] * w1[k1] * w2[k2]

  FREE: w1=w2=gf  (HWE prior [(1-f)^2, 2f(1-f), f^2])
  CHIP: w1=genoProb (from chip genotype), w2=gf
        genoProb: missing->gf; homref->[1-ge,ge/2,ge/2]; het->[ge/2,1-ge,ge/2];
        homalt->[ge/2,ge/2,1-ge];  ge = genoError (default 1e-3)

smLLK = sum_i log(perMarker_i);  f(fMix) = -smLLK.  fMix = 1 - alpha.
FREELK0/CHIPLK0 = f(1.0) (alpha 0); MIX = argmin_alpha f(1-alpha); LK1 = f at min.
Empty markers (no reads) -> perMarker = (sum w1)(sum w2) = 1 -> log 0 (inert).
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pyarrow.parquet as pq
from scipy.optimize import minimize_scalar

PSN_REF = np.array([1.0, 0.5, 0.0])
PSN_ALT = np.array([0.0, 0.5, 1.0])
GENO_ERROR = 1.0e-3


def load(markers_path, pileup_path, max_q=40):
    m = pq.read_table(markers_path)
    ref = np.array(m.column("ref").to_pylist())
    alt = np.array(m.column("alt").to_pylist())
    af = np.asarray(m.column("af"), dtype=np.float64)

    p = pq.read_table(pileup_path)
    marker = np.asarray(p.column("marker"), dtype=np.int64)        # sorted ascending
    base = np.array(p.column("base").to_pylist())
    qual = np.asarray(p.column("qual"), dtype=np.float64)

    e = np.power(10.0, -np.minimum(qual, max_q) / 10.0)
    matchp = 1.0 - e
    other = e / 3.0
    pb_ref = np.where(base == ref[marker], matchp, other)          # P(b | ref allele)
    pb_alt = np.where(base == alt[marker], matchp, other)          # P(b | alt allele)

    uniq, starts = np.unique(marker, return_index=True)            # markers with reads
    afc = np.clip(af[uniq], 0.001, 0.999)
    gf = np.stack([(1 - afc) ** 2, 2 * afc * (1 - afc), afc ** 2], axis=1)  # (nu,3) HWE
    gfo = (gf[:, :, None] * gf[:, None, :]).reshape(len(uniq), 9)            # gf[k1]*gf[k2]
    return dict(n_markers=m.num_rows, pb_ref=pb_ref, pb_alt=pb_alt,
                starts=starts, uniq=uniq, gf=gf, gfo=gfo)


def chip_weights(d, self_geno, geno_error=GENO_ERROR):
    """w1[k1]*w2[k2] for the CHIP model: genoProb(chip) outer HWE, per reads-marker.

    self_geno is a full-panel int8 array (0 missing, 1 homref, 2 het, 3 homalt).
    """
    g = np.asarray(self_geno)[d["uniq"]]
    gf = d["gf"]
    nu = len(g)
    ge = geno_error
    genoprob = np.empty((nu, 3))
    genoprob[g == 1] = [1 - ge, ge / 2, ge / 2]   # homref
    genoprob[g == 2] = [ge / 2, 1 - ge, ge / 2]   # het
    genoprob[g == 3] = [ge / 2, ge / 2, 1 - ge]   # homalt
    miss = g == 0                                  # missing chip geno -> fall back to HWE
    if miss.any():
        genoprob[miss] = gf[miss]
    return (genoprob[:, :, None] * gf[:, None, :]).reshape(nu, 9)  # w1[k1]*w2[k2]


def _seg(fmix, d):
    """markerLK in log space, (nu, 9) -- shared by FREE and CHIP."""
    A = (fmix * PSN_REF[:, None] + (1 - fmix) * PSN_REF[None, :]).ravel()
    B = (fmix * PSN_ALT[:, None] + (1 - fmix) * PSN_ALT[None, :]).ravel()
    baselk = d["pb_ref"][:, None] * A[None, :] + d["pb_alt"][:, None] * B[None, :]
    return np.add.reduceat(np.log(baselk), d["starts"], axis=0)


def neg_llk(fmix, d, w):
    seg = _seg(fmix, d)
    per_marker = (np.exp(seg) * w).sum(axis=1)
    return -np.log(per_marker).sum()


def optimize(d, w, grid=0.05, max_alpha=0.5):
    """Grid scan alpha in [0, max_alpha] of f(1-alpha) = neg_llk, then Brent refine.

    FREE uses max_alpha=0.5 (the model is alpha<->1-alpha symmetric, so the minor
    fraction is the meaningful basin). CHIP fixes the self genotype, breaking that
    symmetry, so it can scan higher and even flag sample swaps (alpha -> 1).
    """
    alphas = np.arange(0.0, max_alpha + 1e-9, grid)
    lks = np.array([neg_llk(1.0 - a, d, w) for a in alphas])
    k = int(lks.argmin())
    llk0 = float(lks[0])  # alpha = 0

    if k == 0:
        lo, hi = 0.0, grid
    elif k == len(alphas) - 1:
        lo, hi = alphas[k - 1], max_alpha
    else:
        lo, hi = alphas[k - 1], alphas[k + 1]
    res = minimize_scalar(lambda a: neg_llk(1.0 - a, d, w), bounds=(lo, hi),
                          method="bounded", options={"xatol": 1e-4})
    alpha_opt, llk1 = float(res.x), float(res.fun)
    if llk1 > llk0:
        alpha_opt, llk1 = 0.0, llk0
    return alpha_opt, llk1, llk0


def _refine_alpha(d, w, alphas, grid_lks):
    """Brent-refine alpha around the grid argmin for one individual's chip weights.
    grid_lks (neg-llk per grid alpha) only locates the bracket; the refined value and
    llk0 are recomputed from w so they're self-consistent. Returns (chipmix,lk1,lk0).
    """
    k = int(grid_lks.argmin())
    lo = alphas[max(k - 1, 0)]
    hi = alphas[min(k + 1, len(alphas) - 1)]
    res = minimize_scalar(lambda a: neg_llk(1.0 - a, d, w), bounds=(lo, hi),
                          method="bounded", options={"xatol": 1e-4})
    llk0 = neg_llk(1.0, d, w)                       # alpha = 0, same w-path
    chipmix, lk1 = float(res.x), float(res.fun)
    if lk1 > llk0:
        chipmix, lk1 = 0.0, llk0
    return chipmix, lk1, llk0


def scan_chip(d, chip_matrix_path, self_id=None, grid=0.05, max_alpha=0.95, top=3):
    """Scan every individual in the chip matrix against the reads (verifyBamID --best,
    Main.cpp:319-353): for each individual fit contamination alpha; the BEST match is
    the lowest-alpha (highest-IBD) individual. Also returns the result for the claimed
    self_id (if genotyped) so callers get both the identity scan AND contamination-vs-
    claimed, plus a ranking of the closest individuals.

    Vectorized: per_marker = sum_k1 genoProb[k1] * G[k1], where G[i,k1] = sum_k2
    exp(seg)[i,k1,k2]*gf[k2] depends only on alpha (not the individual). So at each
    grid alpha we build four per-marker templates (genotype 0/1/2/3) once and gather
    them by the (markers x samples) genotype matrix -- the whole 278-sample x 290k-
    marker scan is then a handful of array ops instead of a Python loop per sample.
    Mathematically identical to optimize(chip_weights(geno)) per individual.

    Returns dict: best=(id,chipmix,lk1,lk0); self=(...)|None; ranking=[(id,alpha,lk1)].
    """
    tbl = pq.read_table(chip_matrix_path)
    sample_ids = tbl.schema.names
    uniq, gf = d["uniq"], d["gf"]                   # reads-markers, HWE prior (nu,3)
    gmat = np.column_stack([np.asarray(tbl.column(s), dtype=np.int8)[uniq]
                            for s in sample_ids])    # (nu, n_samples) genotype codes
    nu, nsamp = gmat.shape
    ge = GENO_ERROR
    gp = np.array([1 - ge, ge / 2, ge / 2,          # rows homref/het/homalt (code 1/2/3)
                   ge / 2, 1 - ge, ge / 2,
                   ge / 2, ge / 2, 1 - ge]).reshape(3, 3)

    alphas = np.arange(0.0, max_alpha + 1e-9, grid)
    lkmat = np.empty((len(alphas), nsamp))
    for ai, a in enumerate(alphas):
        segE = np.exp(_seg(1.0 - a, d)).reshape(nu, 3, 3)
        G = (segE * gf[:, None, :]).sum(axis=2)     # (nu,3): sum_k2 segE*gf
        T = np.empty((nu, 4))
        T[:, 0] = (gf * G).sum(axis=1)              # code 0 (missing) -> HWE on k1
        T[:, 1:] = G @ gp.T                          # codes 1/2/3 -> genoProb on k1
        P = T[np.arange(nu)[:, None], gmat]          # (nu,nsamp) gather by genotype
        lkmat[ai] = -np.log(P).sum(axis=0)
    grid_alpha = lkmat.argmin(axis=0)                # per-sample grid argmin index
    min_alpha = alphas[grid_alpha]
    min_lk = lkmat[grid_alpha, np.arange(nsamp)]
    order = np.lexsort((min_lk, min_alpha))          # best = lowest alpha, then lowest lk

    def _result(j):
        sid = sample_ids[j]
        w = chip_weights(d, np.asarray(tbl.column(sid), dtype=np.int8))
        return (sid, *_refine_alpha(d, w, alphas, lkmat[:, j]))

    best = _result(int(order[0]))
    self_res = _result(sample_ids.index(self_id)) if self_id in sample_ids else None
    ranking = [(sample_ids[j], round(float(min_alpha[j]), 4), round(float(min_lk[j]), 2))
               for j in order[:top]]
    return dict(best=best, self=self_res, ranking=ranking)


def extract_self_geno(chip_vcf, sample_id, markers_path):
    """Per-marker genotype (0 missing/1 homref/2 het/3 homalt) for sample_id, aligned
    to the panel markers, matching verifyBamID setGenotype semantics (half-call's
    missing allele -> REF via atoi('.')==0; first-allele-missing -> whole geno missing).
    """
    from cyvcf2 import VCF
    m = pq.read_table(markers_path)
    chrom = m.column("chrom").to_pylist()
    pos = m.column("pos").to_pylist()
    key_to_idx = {(c, p): i for i, (c, p) in enumerate(zip(chrom, pos))}

    vcf = VCF(chip_vcf, samples=[sample_id])
    if sample_id not in list(vcf.samples):
        # No genotypes for this sample in the chip source -> CHIPMIX is NA by
        # definition (it needs the individual's own genotypes). Caller falls back
        # to FREEMIX. Degrade gracefully rather than aborting a batch run.
        print(f"[chip] sample {sample_id} not in {chip_vcf}; CHIPMIX=NA", file=sys.stderr)
        return None

    geno = np.zeros(m.num_rows, dtype=np.int8)  # default 0 = missing
    for v in vcf:
        idx = key_to_idx.get((v.CHROM, v.POS))
        if idx is None:
            continue
        gt = v.genotype.array()[0]
        a, b = int(gt[0]), int(gt[1])
        if a < 0:                       # first allele missing -> whole genotype missing
            continue
        b = 0 if b < 0 else b           # trailing-missing half-call -> REF
        alt_count = (a > 0) + (b > 0)
        geno[idx] = 1 + alt_count       # homref=1, het=2, homalt=3
    return geno


def main(argv=None):
    ap = argparse.ArgumentParser(description="Estimate FREEMIX/FREELK (+CHIPMIX) from pileup.")
    ap.add_argument("--markers", required=True, help="<prefix>.markers.parquet")
    ap.add_argument("--pileup", required=True, help="<prefix>.pileup.parquet")
    ap.add_argument("--seq-id", default="SAMPLE")
    ap.add_argument("--chip-vcf", help="VCF with the sample's genotypes -> enables CHIPMIX")
    ap.add_argument("--chip-matrix", help="precomputed chip matrix (build-chip) -> CHIPMIX")
    ap.add_argument("--chip-id", help="sample id in chip source (default: --seq-id)")
    ap.add_argument("--best", action="store_true",
                    help="(optional) scan all chip-matrix individuals for the best-matching "
                         "ID -> identity / sample-swap check; requires --chip-matrix")
    ap.add_argument("--max-q", type=int, default=40)
    ap.add_argument("--grid", type=float, default=0.05)
    a = ap.parse_args(argv)

    d = load(a.markers, a.pileup, max_q=a.max_q)
    freemix, freelk1, freelk0 = optimize(d, d["gfo"], grid=a.grid, max_alpha=0.5)

    chip = ["NA", "NA", "NA"]
    chip_id_out = "NA"
    if a.best:
        if not a.chip_matrix:
            sys.exit("--best requires --chip-matrix (it scans all individuals' genotypes)")
        self_id = a.chip_id or a.seq_id
        res = scan_chip(d, a.chip_matrix, self_id=self_id, grid=a.grid)
        best_id = res["best"][0]
        if res["self"] is not None:               # CHIPMIX columns = vs CLAIMED self
            _, cmix, clk1, clk0 = res["self"]
            chip = [f"{cmix:.5f}", f"{clk1:.2f}", f"{clk0:.2f}"]
            chip_id_out = self_id
            swap = "NO" if best_id == self_id else "YES"
        else:
            swap = "NA"                            # claimed id not genotyped -> unconfirmable
        rank = "; ".join(f"{i}:{al}" for i, al, _ in res["ranking"])
        print(f"BEST_MATCH={best_id}  SWAP={swap}  best_chipmix={res['best'][1]:.5f}  "
              f"self_chipmix={chip[0]}  top={rank}", file=sys.stderr)
    elif a.chip_matrix or a.chip_vcf:
        if a.chip_matrix:
            from . import build_chip
            sg = build_chip.load_sample(a.chip_matrix, a.chip_id or a.seq_id)
        else:
            sg = extract_self_geno(a.chip_vcf, a.chip_id or a.seq_id, a.markers)
        if sg is not None:
            w = chip_weights(d, sg)
            chipmix, chiplk1, chiplk0 = optimize(d, w, grid=a.grid, max_alpha=0.95)
            chip = [f"{chipmix:.5f}", f"{chiplk1:.2f}", f"{chiplk0:.2f}"]
            chip_id_out = a.chip_id or a.seq_id
            print(f"CHIPMIX={chip[0]}  CHIPLK1={chip[1]}  CHIPLK0={chip[2]}", file=sys.stderr)

    print(f"FREEMIX={freemix:.5f}  FREELK1={freelk1:.2f}  FREELK0={freelk0:.2f}",
          file=sys.stderr)
    cols = [a.seq_id, "ALL", chip_id_out, str(d["n_markers"]), "NA", "NA",
            f"{freemix:.5f}", f"{freelk1:.2f}", f"{freelk0:.2f}", "NA", "NA",
            chip[0], chip[1], chip[2], "NA", "NA", "NA", "NA", "NA"]
    print("\t".join(cols))


if __name__ == "__main__":
    main()
