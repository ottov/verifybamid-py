"""Load a VCPA marker table as estimator input (estimate --table).

The VCPA marker pileup writes one row per covered panel marker:

    chrom  pos  ref  alt  refDepth  altDepth  refQsum  altQsum

where Qsum is the phred sum over that allele's reads. Each allele becomes one
weighted observation at its mean quality round(Qsum/depth) -- the same rule
ad2pileup_q.py used to expand tables into per-read pileups, so estimates match the
expanded path while doing ~20x less work.

A source is a full table (*.ad.tsv.gz), or a shard set as marker-shards.sh writes
it: a <SM>.markers.manifest.tsv (or the directory holding exactly one) listing each
shard's file, panel_markers, rows and bytes; or the same on S3, as an
s3://bucket/.../markers/ prefix or the manifest's key. A shard set is read only
through its manifest and must match it exactly -- a missing or altered shard would
otherwise silently drop part of a chromosome from the estimate.
"""

from __future__ import annotations

import gzip
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.parquet as pq

from . import estimate

COLUMNS = ["chrom", "pos", "ref", "alt", "rd", "ad", "rq", "aq"]
_TYPES = dict(chrom=pa.string(), pos=pa.int64(), ref=pa.string(), alt=pa.string(),
              rd=pa.int64(), ad=pa.int64(), rq=pa.int64(), aq=pa.int64())
MANIFEST_SUFFIX = ".markers.manifest.tsv"


class TableError(Exception):
    """The source is not a complete, consistent marker table for this panel."""


def parse(text: bytes) -> pa.Table:
    """Table rows from uncompressed TSV bytes."""
    if not text:
        return pa.table({c: pa.array([], t) for c, t in _TYPES.items()})
    return pcsv.read_csv(
        pa.BufferReader(text),
        read_options=pcsv.ReadOptions(column_names=COLUMNS),
        parse_options=pcsv.ParseOptions(delimiter="\t"),
        convert_options=pcsv.ConvertOptions(column_types=_TYPES))


def _marker_index(m, chrom, pos):
    """Panel index of each (chrom, pos), or -1 where it is not a panel marker."""
    names = {c: k for k, c in enumerate(dict.fromkeys(m.column("chrom").to_pylist()))}
    pkey = (np.array([names[c] for c in m.column("chrom").to_pylist()], dtype=np.int64) << 32) \
        | np.asarray(m.column("pos"), dtype=np.int64)
    order = np.argsort(pkey, kind="stable")
    code = np.array([names.get(c, -1) for c in chrom], dtype=np.int64)
    key = (code << 32) | pos
    at = np.searchsorted(pkey[order], key)
    at = np.minimum(at, len(order) - 1)
    hit = (code >= 0) & (pkey[order][at] == key)
    return np.where(hit, order[at], -1)


def observations(m, rows: pa.Table):
    """(marker, base, qual, weight) sorted by marker -- one entry per non-zero allele --
    and the number of rows that are not panel markers (ignored)."""
    idx = _marker_index(m, rows.column("chrom").to_pylist(), np.asarray(rows.column("pos")))
    keep = idx >= 0
    idx = idx[keep]
    col = {c: np.asarray(rows.column(c))[keep] for c in ("rd", "ad", "rq", "aq")}
    ref = np.array(rows.column("ref").to_pylist())[keep]
    alt = np.array(rows.column("alt").to_pylist())[keep]
    has_r, has_a = col["rd"] > 0, col["ad"] > 0
    marker = np.concatenate([idx[has_r], idx[has_a]])
    base = np.concatenate([ref[has_r], alt[has_a]])
    qual = np.concatenate([np.round(col["rq"][has_r] / col["rd"][has_r]),
                           np.round(col["aq"][has_a] / col["ad"][has_a])])
    weight = np.concatenate([col["rd"][has_r], col["ad"][has_a]]).astype(np.float64)
    o = np.argsort(marker, kind="stable")
    return marker[o], base[o], qual[o], weight[o], int((~keep).sum())


def _gunzip(data: bytes, where: str) -> bytes:
    try:
        return gzip.decompress(data)
    except (OSError, EOFError) as e:                 # not gzip, or cut mid-stream
        raise TableError(f"{where}: {e}") from None


def read_manifest(text: str, where: str):
    """(declared panel marker count, [shard entries]) from manifest text."""
    head = [ln for ln in text.splitlines() if ln.startswith("#")]
    declared = next((int(w.split("=", 1)[1]) for ln in head for w in ln[1:].split()
                     if w.startswith("markers=")), None)
    body = [ln.split("\t") for ln in text.splitlines() if ln and not ln.startswith("#")]
    if declared is None or not body or body[0][:2] != ["shard", "file"]:
        raise TableError(f"{where}: not a marker-shards manifest")
    cols = body[0]
    shards = [dict(zip(cols, r)) for r in body[1:]]
    for e in shards:
        for k in ("panel_markers", "rows", "bytes"):
            e[k] = int(e[k])
    return declared, shards


def join_shards(manifest_text: str, where: str, fetch, n_panel: int) -> bytes:
    """Concatenated table text of a shard set, after checking every shard against the
    manifest and the manifest against the panel. fetch(name) returns a shard's bytes
    or raises FileNotFoundError."""
    declared, shards = read_manifest(manifest_text, where)
    if declared != n_panel:
        raise TableError(f"{where}: manifest is for a {declared}-marker panel, "
                         f"markers file has {n_panel}")
    listed = sum(e["panel_markers"] for e in shards)
    if listed != declared:
        raise TableError(f"{where}: shards cover {listed} panel markers, manifest "
                         f"declares {declared} -- manifest is incomplete")
    parts = []
    for e in shards:
        try:
            data = fetch(e["file"])
        except FileNotFoundError:
            raise TableError(f"{where}: shard {e['file']} is missing") from None
        if len(data) != e["bytes"]:
            raise TableError(f"{where}: shard {e['file']} has {len(data)} bytes, "
                             f"manifest says {e['bytes']}")
        text = _gunzip(data, f"{where}: shard {e['file']}")
        n_rows = text.count(b"\n")
        if n_rows != e["rows"]:
            raise TableError(f"{where}: shard {e['file']} has {n_rows} rows, "
                             f"manifest says {e['rows']}")
        parts.append(text)
    return b"".join(parts)


def _local_manifest(path: Path) -> Path | None:
    if path.is_dir():
        found = sorted(path.glob("*" + MANIFEST_SUFFIX))
        if len(found) != 1:
            raise TableError(f"{path}: expected exactly one *{MANIFEST_SUFFIX}, "
                             f"found {len(found)}")
        return found[0]
    return path if path.name.endswith(MANIFEST_SUFFIX) else None


class S3Store:
    """list/get over boto3, with botocore's standard retry mode. A missing key raises
    FileNotFoundError, as a missing local shard does."""

    def __init__(self, client=None):
        if client is None:
            import boto3
            from botocore.config import Config
            client = boto3.client("s3", config=Config(retries={"max_attempts": 10,
                                                               "mode": "standard"}))
        self.client = client

    def list(self, bucket, prefix):
        pages = self.client.get_paginator("list_objects_v2").paginate(Bucket=bucket, Prefix=prefix)
        return [o["Key"] for page in pages for o in page.get("Contents", [])]

    def get(self, bucket, key):
        try:
            return self.client.get_object(Bucket=bucket, Key=key)["Body"].read()
        except self.client.exceptions.NoSuchKey:
            raise FileNotFoundError(f"s3://{bucket}/{key}") from None


def _read_s3(uri: str, n_panel: int, s3) -> bytes:
    bucket, _, key = uri[len("s3://"):].partition("/")
    if key.endswith(MANIFEST_SUFFIX):
        manifest_key = key
    else:
        prefix = key.rstrip("/") + "/"
        found = [k for k in s3.list(bucket, prefix)
                 if k.endswith(MANIFEST_SUFFIX) and "/" not in k[len(prefix):]]
        if len(found) != 1:
            raise TableError(f"{uri}: expected exactly one *{MANIFEST_SUFFIX}, found {len(found)}")
        manifest_key = found[0]
    where = f"s3://{bucket}/{manifest_key}"
    try:
        text = s3.get(bucket, manifest_key).decode()
    except FileNotFoundError:
        raise TableError(f"{where}: manifest not found") from None
    base = manifest_key[:manifest_key.rfind("/") + 1]
    names = [e["file"] for e in read_manifest(text, where)[1]]
    with ThreadPoolExecutor(max_workers=16) as pool:
        futures = {n: pool.submit(s3.get, bucket, base + n) for n in names}
    fetched = {}
    for n, f in futures.items():
        try:
            fetched[n] = f.result()
        except FileNotFoundError:
            pass                                  # join_shards reports it by name

    def fetch(name):
        if name not in fetched:
            raise FileNotFoundError(name)
        return fetched[name]
    return join_shards(text, where, fetch, n_panel)


def read_text(source, n_panel: int, s3=None) -> bytes:
    """Uncompressed table text for a full table, manifest, shard directory or S3 shard set."""
    if str(source).startswith("s3://"):
        return _read_s3(str(source), n_panel, s3 or S3Store())
    path = Path(source)
    manifest = _local_manifest(path)
    if manifest is None:
        try:
            data = path.read_bytes()
        except OSError as e:
            raise TableError(f"{path}: {e.strerror}") from None
        return _gunzip(data, str(path))
    return join_shards(manifest.read_text(), str(manifest),
                       lambda name: (manifest.parent / name).read_bytes(), n_panel)


def load(markers_path, source, max_q=40, s3=None):
    m = pq.read_table(markers_path)
    rows = parse(read_text(source, m.num_rows, s3=s3))
    marker, base, qual, weight, stray = observations(m, rows)
    d = estimate.prepare(m, marker, base, qual, weight=weight, max_q=max_q)
    d.update(n_reads=int(weight.sum()), covered=len(d["uniq"]), not_in_panel=stray)
    return d
