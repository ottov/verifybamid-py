"""Stage 1: stream a CRAM (local or S3) and pile up bases at panel markers.

Replicates verifyBamID 1.1.3 BamPileBases::readMarker (git/.../BamPileBases.cpp:121-227)
exactly, so #READS / AVG_DP / the depth histogram match the original .selfSM/.depthSM.

Engine: a single sequential pass per marker chunk (one fetch over the chunk's span),
tallying each read into the markers it overlaps via a CIGAR walk. This decodes each
CRAM slice once and issues large sequential range reads -- essential for S3 streaming,
where per-marker random fetches would mean one network round-trip per marker. Reads
arrive in coordinate (= file) order, which is exactly the order verifyBamID uses, so
the "first mate claims the read name" dedup is reproduced faithfully.

Per-read filters (in order), matching readMarker:
  * mapping quality >= minMapQ (default 10)
  * (flag & 0x0704) == 0  = exclude unmapped|secondary|QCfail|duplicate
      (supplementary 0x800 is NOT excluded, matching verifyBamID)
  * the read name is claimed into the per-marker dedup set the moment the read passes
    mapQ+flags -- BEFORE the base passes (BamPileBases.cpp:208); the other mate is dropped
  * marker must align to a query base (CIGAR M/=/X); a deletion / ref-skip claims the
    name but yields no base
  * base != 'N' and base quality (phred) >= minQ (default 13)
  * at most maxDepth (default 20) bases per marker
maxQ (40) is NOT applied here; it caps quality in the likelihood (Stage 2).
"""

from __future__ import annotations

import argparse
import bisect
import os
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pysam

warnings.filterwarnings("ignore", message="multiple_iterators")

EXCLUDE_FLAG = 0x0704  # unmapped | secondary | qcfail | duplicate
DEFAULTS = dict(min_mapq=10, min_q=13, max_q=40, max_depth=20)


def open_cram(cram_url: str, ref: str, crai_url: str | None):
    if cram_url.startswith(("http://", "https://", "s3://")):
        os.environ.setdefault("CURL_CA_BUNDLE", "/etc/ssl/certs/ca-certificates.crt")
    kw = {"index_filename": crai_url} if crai_url else {}
    return pysam.AlignmentFile(cram_url, "rc", reference_filename=ref, **kw)


def base_at(cigartuples, seq, quals, rstart, refpos0):
    """Query (base, qual) aligned to 0-based reference position refpos0, or None if
    the position is a deletion/ref-skip in this read (CigarRoller INDEX_NA)."""
    pos = rstart  # 0-based reference cursor
    qi = 0        # query index into full seq (incl. soft-clips)
    for op, ln in cigartuples:
        if op == 0 or op == 7 or op == 8:        # M/=/X: consume ref + query
            if pos <= refpos0 < pos + ln:
                o = qi + (refpos0 - pos)
                return seq[o], quals[o]
            pos += ln
            qi += ln
        elif op == 1 or op == 4:                 # I/S: consume query only
            qi += ln
        elif op == 2 or op == 3:                 # D/N: consume ref only
            if pos <= refpos0 < pos + ln:
                return None
            pos += ln
        # H(5)/P(6): consume nothing
    return None


def pile_span(af, contig, positions, gidx, *, min_mapq, min_q, max_q, max_depth):
    """Pile a contiguous run of markers on one contig in a single sequential pass."""
    mpos0 = [p - 1 for p in positions]  # 0-based, ascending
    n = len(positions)
    seen = [None] * n
    bases = [None] * n
    quals = [None] * n
    full = bytearray(n)

    for read in af.fetch(contig, mpos0[0], positions[-1]):
        if read.mapping_quality < min_mapq:
            continue
        if read.flag & EXCLUDE_FLAG:
            continue
        re = read.reference_end
        if re is None:
            continue
        rs = read.reference_start
        lo = bisect.bisect_left(mpos0, rs)        # markers with mpos0 >= rs
        hi = bisect.bisect_left(mpos0, re)        # markers with mpos0 <  re
        if lo == hi:
            continue
        name = read.query_name
        cig = read.cigartuples
        seq = read.query_sequence
        rq = read.query_qualities
        for k in range(lo, hi):
            if full[k]:
                continue
            s = seen[k]
            if s is None:
                s = seen[k] = set()
            if name in s:                          # mate already claimed this site
                continue
            s.add(name)                            # claim before the base is tested
            bq = base_at(cig, seq, rq, rs, mpos0[k])
            if bq is None:                         # deletion/ref-skip: claimed, no base
                continue
            b, q = bq
            if q < min_q or b == "N":
                continue
            bl = bases[k]
            if bl is None:
                bl = bases[k] = []
                quals[k] = []
            bl.append(b)
            quals[k].append(q)
            if len(bl) >= max_depth:
                full[k] = 1

    out_g, out_d, mi, bb, qq = [], [], [], [], []
    for k in range(n):
        d = 0 if bases[k] is None else len(bases[k])
        out_g.append(gidx[k])
        out_d.append(d)
        if d:
            mi.extend([gidx[k]] * d)
            bb.extend(bases[k])
            qq.extend(quals[k])
    return out_g, out_d, mi, bb, qq


# ---- worker (one chunk = contiguous markers on one contig) -----------------
_W = {}


def _init(cram, ref, crai, filt):
    _W.update(cram=cram, ref=ref, crai=crai, filt=filt)
    _W["af"] = open_cram(cram, ref, crai)


def _work(chunk, retries=5):
    contig, positions, gidx = chunk
    last = None
    for attempt in range(retries):
        try:
            return pile_span(_W["af"], contig, positions, gidx, **_W["filt"])
        except (OSError, ValueError) as e:
            # transient remote read (truncated range, reset connection, etc.):
            # reopen the handle and retry. Re-piling a chunk is idempotent.
            last = e
            try:
                _W["af"].close()
            except Exception:
                pass
            time.sleep(min(2 ** attempt, 10))
            _W["af"] = open_cram(_W["cram"], _W["ref"], _W["crai"])
    raise RuntimeError(f"chunk {contig}:{positions[0]}-{positions[-1]} failed "
                       f"after {retries} retries: {last}")


def run(cram_url, ref, panel_path, *, crai_url=None, out=None, contigs=None,
        jobs=1, max_span=None, **filt):
    f = {**DEFAULTS, **{k: v for k, v in filt.items() if v is not None}}
    panel = pq.read_table(panel_path)
    n_snps = panel.num_rows  # ALL markers count toward #SNPS / AVG_DP denominator
    chrom = panel.column("chrom").to_pylist()
    pos = panel.column("pos").to_pylist()

    af0 = open_cram(cram_url, ref, crai_url)
    cram_contigs = set(af0.references)
    af0.close()
    want = set(contigs) if contigs else None

    # group markers by contig (panel is position-sorted within contig)
    by_contig: dict[str, tuple[list, list]] = {}
    for i in range(n_snps):
        c = chrom[i]
        if c in cram_contigs and (want is None or c in want):
            ps, gs = by_contig.setdefault(c, ([], []))
            ps.append(pos[i])
            gs.append(i)
    n_work = sum(len(v[0]) for v in by_contig.values())

    # contiguous per-contig chunks, sized so each job gets ~8 chunks. max_span caps
    # the genomic span per chunk: leave it None for a dense panel (one big sequential
    # read), set it (~100-200kb) for a sparse/downsampled panel so each fetch pulls
    # only the slice(s) around its markers instead of everything in between.
    chunk_target = max(1000, n_work // max(jobs * 8, 1))
    chunks = []
    for c, (ps, gs) in by_contig.items():
        i = 0
        while i < len(ps):
            j = i + 1
            while (j < len(ps) and (j - i) < chunk_target
                   and (max_span is None or ps[j] - ps[i] <= max_span)):
                j += 1
            chunks.append((c, ps[i:j], gs[i:j]))
            i = j
    print(f"piling {n_work:,} markers in {len(chunks)} chunks on {jobs} job(s) "
          f"(of {n_snps:,} #SNPS)", file=sys.stderr)

    depth = np.zeros(n_snps, dtype=np.int32)
    rows_mi, rows_b, rows_q = [], [], []
    t0 = time.time()

    def absorb(res):
        g, d, mi, b, q = res
        for gi, di in zip(g, d):
            depth[gi] = di
        rows_mi.extend(mi); rows_b.extend(b); rows_q.extend(q)

    if jobs == 1:
        _init(cram_url, ref, crai_url, f)
        for n, ch in enumerate(chunks, 1):
            absorb(_work(ch))
            if n % 10 == 0 or n == len(chunks):
                print(f"  {n}/{len(chunks)} chunks ({time.time()-t0:.0f}s)", file=sys.stderr)
    else:
        with Pool(jobs, initializer=_init, initargs=(cram_url, ref, crai_url, f)) as pool:
            for n, res in enumerate(pool.imap_unordered(_work, chunks), 1):
                absorb(res)
                if n % 10 == 0 or n == len(chunks):
                    print(f"  {n}/{len(chunks)} chunks ({time.time()-t0:.0f}s)", file=sys.stderr)

    total_reads = int(depth.sum())
    avg_dp = total_reads / n_snps if n_snps else 0.0
    hist = np.bincount(np.minimum(depth, f["max_depth"]), minlength=f["max_depth"] + 1)

    print(f"\n#SNPS={n_snps:,}  #READS={total_reads:,}  AVG_DP={avg_dp:.2f}  "
          f"({time.time()-t0:.0f}s)", file=sys.stderr)
    print("depth histogram (depth: nSNPs):", file=sys.stderr)
    for d in range(f["max_depth"], -1, -1):
        if hist[d]:
            print(f"  {d:>2}: {hist[d]:,}", file=sys.stderr)

    if out:
        order = np.argsort(rows_mi, kind="stable") if rows_mi else []
        pq.write_table(
            pa.table({"marker": pa.array(np.asarray(rows_mi)[order], pa.int32()),
                      "base": pa.array(np.asarray(rows_b)[order], pa.string()),
                      "qual": pa.array(np.asarray(rows_q)[order], pa.int8())}),
            f"{out}.pileup.parquet")
        pq.write_table(panel.append_column("depth", pa.array(depth, pa.int32())),
                       f"{out}.markers.parquet")
        print(f"\nwrote {out}.pileup.parquet and {out}.markers.parquet", file=sys.stderr)

    return dict(n_snps=n_snps, n_reads=total_reads, avg_dp=avg_dp)


def main(argv=None):
    p = argparse.ArgumentParser(description="Stream CRAM and pile up at panel markers.")
    p.add_argument("--cram", required=True, help="local path or presigned http(s) URL")
    p.add_argument("--crai", help="presigned URL / path for the .crai")
    p.add_argument("--ref", required=True)
    p.add_argument("--panel", required=True, help="panel.parquet from build_panel")
    p.add_argument("--out", help="output prefix (.pileup.parquet, .markers.parquet)")
    p.add_argument("--contig", action="append", dest="contigs",
                   help="restrict to contig(s); repeatable")
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--max-span", type=int,
                   help="cap genomic span (bp) per fetch; set ~150000 for sparse panels")
    p.add_argument("--min-mapq", type=int)
    p.add_argument("--min-q", type=int)
    p.add_argument("--max-q", type=int)
    p.add_argument("--max-depth", type=int)
    a = p.parse_args(argv)
    run(a.cram, a.ref, a.panel, crai_url=a.crai, out=a.out, contigs=a.contigs,
        jobs=a.jobs, max_span=a.max_span, min_mapq=a.min_mapq, min_q=a.min_q,
        max_q=a.max_q, max_depth=a.max_depth)


if __name__ == "__main__":
    main()
