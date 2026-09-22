"""Synthetic panels, marker tables and shard sets for the tests.

Real marker tables are individual-level genotype data and never go in the repo, so
the tests build small ones here in the same formats: the 8-column table written by
the VCPA marker pileup, the shard set + manifest written by marker-shards.sh, and
the per-read pileup that ad2pileup_q.py expanded tables into before --table existed.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

BASES = np.array(list("ACGT"))


def make_panel(path: Path, n_per_chrom=(180, 120), seed=1):
    """markers.parquet with chrom/pos/ref/alt/af; returns the columns as lists."""
    rng = np.random.default_rng(seed)
    chrom, pos, ref, alt = [], [], [], []
    for c, n in enumerate(n_per_chrom, start=1):
        p = np.sort(rng.choice(np.arange(10_000, 50_000_000), n, replace=False))
        for x in p:
            r, a = rng.choice(BASES, 2, replace=False)
            chrom.append(f"chr{c}"); pos.append(int(x)); ref.append(str(r)); alt.append(str(a))
    af = rng.uniform(0.05, 0.95, len(pos)).round(4).tolist()
    pq.write_table(pa.table({"chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "af": af}), path)
    return dict(chrom=chrom, pos=pos, ref=ref, alt=alt, af=af)


def make_rows(panel, alpha=0.08, depth=30, skip_every=0, seed=2):
    """Table rows (chrom, pos, ref, alt, rd, ad, rq, aq) for a sample contaminated at
    alpha. Qualities vary per read so Qsum/depth is usually not a whole number.
    skip_every=k drops every k-th marker (uncovered), as real tables omit them."""
    rng = np.random.default_rng(seed)
    rows = []
    for i, (c, p, r, a, f) in enumerate(zip(panel["chrom"], panel["pos"], panel["ref"],
                                            panel["alt"], panel["af"])):
        if skip_every and i % skip_every == 0:
            continue
        g_self, g_cont = rng.binomial(2, f), rng.binomial(2, f)
        n = max(1, rng.poisson(depth))
        from_cont = rng.random(n) < alpha
        p_alt = np.where(from_cont, g_cont / 2, g_self / 2)
        is_alt = rng.random(n) < p_alt
        q = rng.integers(18, 42, n)
        rd, ad = int((~is_alt).sum()), int(is_alt.sum())
        rows.append((c, p, r, a, rd, ad, int(q[~is_alt].sum()), int(q[is_alt].sum())))
    return rows


def table_text(rows) -> bytes:
    return "".join("\t".join(map(str, r)) + "\n" for r in rows).encode()


def write_table(path: Path, rows):
    path.write_bytes(gzip.compress(table_text(rows)))
    return path


def write_expanded_pileup(path: Path, panel, rows):
    """What ad2pileup_q.py did: one pileup row per read, every read of an allele at
    that allele's mean quality round(Qsum/depth)."""
    idx = {(c, p): i for i, (c, p) in enumerate(zip(panel["chrom"], panel["pos"]))}
    mi, bs, qs = [], [], []
    for c, p, r, a, rd, ad, rq, aq in rows:
        i = idx.get((c, p))
        if i is None:
            continue
        if rd:
            mi += [i] * rd; bs += [r] * rd; qs += [int(round(rq / rd))] * rd
        if ad:
            mi += [i] * ad; bs += [a] * ad; qs += [int(round(aq / ad))] * ad
    o = np.argsort(np.asarray(mi, dtype=np.int64), kind="stable")
    pq.write_table(pa.table({
        "marker": pa.array(np.asarray(mi, dtype=np.int32)[o], pa.int32()),
        "base": pa.array(np.asarray(bs)[o].tolist(), pa.string()),
        "qual": pa.array(np.asarray(qs, dtype=np.int8)[o], pa.int8()),
    }), path)
    return path


def shard_files(panel, rows, sm="SM", cap=50):
    """{file name: gz bytes} plus manifest text, laid out as marker-shards.sh does:
    each chromosome cut into the fewest near-equal runs of panel markers under cap,
    table rows routed by position, every shard written even if empty."""
    by_chrom = {}
    for i, c in enumerate(panel["chrom"]):
        by_chrom.setdefault(c, []).append(i)
    shards = []                                   # (label, [panel indices])
    for c, ix in by_chrom.items():
        k = -(-len(ix) // cap)
        per = len(ix) / k
        for j in range(k):
            part = [x for n, x in enumerate(ix) if int(n / per) == j]
            shards.append((f"{c}_{j + 1}of{k}", part))
    row_at = {(r[0], r[1]): r for r in rows}
    files, lines = {}, []
    for s, (lab, part) in enumerate(shards, start=1):
        name = f"{sm}.markers.s{s:03d}.{lab}.ad.tsv.gz"
        mine = [row_at[(panel["chrom"][i], panel["pos"][i])] for i in part
                if (panel["chrom"][i], panel["pos"][i]) in row_at]
        files[name] = gzip.compress(table_text(mine))
        lines.append([s, name, panel["chrom"][part[0]], panel["pos"][part[0]],
                      panel["pos"][part[-1]], len(part), len(mine), len(files[name])])
    manifest = (f"# sample={sm} markers={len(panel['pos'])} max_markers_per_shard={cap} "
                f"panel=synthetic md5:000000000000 min_maf:0.05\n"
                "# rows: chrom pos ref alt refDepth altDepth refQsum altQsum. "
                "zcat the shards in name order for the full table.\n"
                "shard\tfile\tchrom\tfirst_pos\tlast_pos\tpanel_markers\trows\tbytes\n" +
                "".join("\t".join(map(str, x)) + "\n" for x in lines))
    return files, manifest


def write_shards(directory: Path, panel, rows, sm="SM", cap=50):
    directory.mkdir(parents=True, exist_ok=True)
    files, manifest = shard_files(panel, rows, sm=sm, cap=cap)
    for name, data in files.items():
        (directory / name).write_bytes(data)
    (directory / f"{sm}.markers.manifest.tsv").write_text(manifest)
    return directory / f"{sm}.markers.manifest.tsv"
