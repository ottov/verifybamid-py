"""--table: estimate straight from a VCPA marker table or its shard set."""

from __future__ import annotations

import pytest

import synth
from verifybamid_py import estimate, table

ALPHAS = (0.0, 0.01, 0.07, 0.2, 0.45)


@pytest.fixture
def sample(tmp_path):
    panel = synth.make_panel(tmp_path / "m.parquet")
    rows = synth.make_rows(panel, alpha=0.08, skip_every=7)
    return dict(markers=tmp_path / "m.parquet", panel=panel, rows=rows, dir=tmp_path)


def selfsm_row(capsys, argv):
    estimate.main(argv)
    return capsys.readouterr().out.rstrip("\n").split("\t")


def test_pileup_path_output_is_unchanged(sample, capsys):
    # Pinned from 5c5a2ba before --table existed: the per-read path must not move.
    p = synth.write_expanded_pileup(sample["dir"] / "p.parquet", sample["panel"], sample["rows"])
    row = selfsm_row(capsys, ["--markers", str(sample["markers"]), "--pileup", str(p),
                              "--seq-id", "S"])
    assert row == ["S", "ALL", "NA", "300", "NA", "NA", "0.07314", "2672.47", "2945.71",
                   "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]


def test_table_likelihood_equals_the_expanded_per_read_likelihood(sample):
    p = synth.write_expanded_pileup(sample["dir"] / "p.parquet", sample["panel"], sample["rows"])
    t = synth.write_table(sample["dir"] / "t.ad.tsv.gz", sample["rows"])
    d_reads = estimate.load(sample["markers"], p)
    d_table = table.load(sample["markers"], t)
    for a in ALPHAS:
        want = estimate.neg_llk(1 - a, d_reads, d_reads["gfo"])
        got = estimate.neg_llk(1 - a, d_table, d_table["gfo"])
        assert got == pytest.approx(want, rel=1e-12), a


def test_rows_off_the_panel_are_counted_and_do_not_change_the_estimate(sample):
    stray = [("chr1", 7, "A", "C", 30, 0, 900, 0), ("chrUn", 5000, "G", "T", 10, 10, 300, 300)]
    clean = table.load(sample["markers"], synth.write_table(sample["dir"] / "a.ad.tsv.gz", sample["rows"]))
    mixed = table.load(sample["markers"],
                       synth.write_table(sample["dir"] / "b.ad.tsv.gz", sample["rows"] + stray))
    assert mixed["not_in_panel"] == 2
    assert estimate.neg_llk(0.9, mixed, mixed["gfo"]) == estimate.neg_llk(0.9, clean, clean["gfo"])


def test_load_reports_reads_and_covered_markers(sample):
    # 300 markers, every 7th dropped -> 257 covered; reads counted from the fixture rows
    d = table.load(sample["markers"], synth.write_table(sample["dir"] / "t.ad.tsv.gz", sample["rows"]))
    assert d["covered"] == 257
    assert d["n_reads"] == 7611
    assert d["n_markers"] == 300


# ---- estimate --table

def test_cli_table_row_matches_the_per_read_estimate_and_fills_reads_and_depth(sample, capsys):
    t = synth.write_table(sample["dir"] / "t.ad.tsv.gz", sample["rows"])
    row = selfsm_row(capsys, ["--markers", str(sample["markers"]), "--table", str(t),
                              "--seq-id", "S"])
    # FREEMIX/LKs as pinned for the per-read path above; 7611 reads / 300 #SNPS = 25.37
    assert row == ["S", "ALL", "NA", "300", "7611", "25.37", "0.07314", "2672.47", "2945.71",
                   "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"]


@pytest.mark.parametrize("n_rows", [0, 100])
def test_cli_refuses_a_table_that_covers_too_little_of_the_panel(sample, capsys, n_rows):
    # 100 of 300 markers is 0.33 coverage, under the default --min-cov-frac 0.80
    t = synth.write_table(sample["dir"] / "t.ad.tsv.gz", sample["rows"][:n_rows])
    with pytest.raises(SystemExit) as e:
        estimate.main(["--markers", str(sample["markers"]), "--table", str(t)])
    assert "covered" in str(e.value.code)
    assert capsys.readouterr().out == ""


def test_cli_min_cov_frac_sets_the_coverage_floor(sample, capsys):
    t = synth.write_table(sample["dir"] / "t.ad.tsv.gz", sample["rows"][:100])
    row = selfsm_row(capsys, ["--markers", str(sample["markers"]), "--table", str(t),
                              "--min-cov-frac", "0.3"])
    assert row[4] == str(sum(r[4] + r[5] for r in sample["rows"][:100]))


def test_cli_reports_a_broken_shard_set_without_a_row(sample, capsys):
    manifest = synth.write_shards(sample["dir"] / "shards", sample["panel"], sample["rows"])
    next(manifest.parent.glob("*.s001.*")).unlink()
    with pytest.raises(SystemExit) as e:
        estimate.main(["--markers", str(sample["markers"]), "--table", str(manifest)])
    assert ".s001." in str(e.value.code)
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("content", [None, b"chr1\t817341\tA\tG\t2\t40\t60\t1230\n"])
def test_cli_reports_an_unreadable_table_without_a_row(sample, capsys, content):
    t = sample["dir"] / "t.ad.tsv.gz"                  # missing, or plain text named .gz
    if content is not None:
        t.write_bytes(content)
    with pytest.raises(SystemExit) as e:
        estimate.main(["--markers", str(sample["markers"]), "--table", str(t)])
    assert str(t) in str(e.value.code)
    assert capsys.readouterr().out == ""


# ---- shard sets (marker-shards.sh layout: <SM>.markers.sNNN.<chrom>_<j>of<k>.ad.tsv.gz + manifest)

def llks(d):
    return [estimate.neg_llk(1 - a, d, d["gfo"]) for a in ALPHAS]


@pytest.fixture
def shards(sample):
    manifest = synth.write_shards(sample["dir"] / "shards", sample["panel"], sample["rows"], cap=50)
    full = table.load(sample["markers"], synth.write_table(sample["dir"] / "t.ad.tsv.gz", sample["rows"]))
    return dict(sample, manifest=manifest, full=full)


def test_manifest_loads_the_same_sample_as_the_full_table(shards):
    d = table.load(shards["markers"], shards["manifest"])
    assert llks(d) == llks(shards["full"])
    assert (d["n_reads"], d["covered"]) == (shards["full"]["n_reads"], shards["full"]["covered"])


def test_directory_holding_one_manifest_is_read_as_its_shard_set(shards):
    d = table.load(shards["markers"], shards["manifest"].parent)
    assert llks(d) == llks(shards["full"])


def test_missing_shard_fails_and_names_it(shards):
    gone = sorted(shards["manifest"].parent.glob("*.s002.*"))[0]
    gone.unlink()
    with pytest.raises(table.TableError, match=gone.name):
        table.load(shards["markers"], shards["manifest"])


def test_shard_that_differs_from_its_manifest_entry_fails(shards):
    victim = sorted(shards["manifest"].parent.glob("*.s003.*"))[0]
    rows = [r for r in shards["rows"] if r[0] == "chr1"][100:101]      # a different, shorter shard
    synth.write_table(victim, rows)
    with pytest.raises(table.TableError, match=victim.name):
        table.load(shards["markers"], shards["manifest"])


def test_shard_with_the_recorded_size_but_other_row_count_fails(shards):
    lines = shards["manifest"].read_text().splitlines(keepends=True)
    f = lines[5].rstrip("\n").split("\t")                           # 3rd shard's entry
    f[6] = str(int(f[6]) + 1)                                       # rows: one more than written
    lines[5] = "\t".join(f) + "\n"
    shards["manifest"].write_text("".join(lines))
    with pytest.raises(table.TableError, match=f"{f[1]} has {int(f[6]) - 1} rows"):
        table.load(shards["markers"], shards["manifest"])


def test_truncated_shard_fails_and_names_it(shards):
    victim = sorted(shards["manifest"].parent.glob("*.s002.*"))[0]
    victim.write_bytes(victim.read_bytes()[:-30])                   # cut mid-stream
    with pytest.raises(table.TableError, match=victim.name):
        table.load(shards["markers"], shards["manifest"])


def test_manifest_missing_an_entry_fails(shards):
    lines = shards["manifest"].read_text().splitlines(keepends=True)
    shards["manifest"].write_text("".join(lines[:-1]))                # drop the last shard's line
    with pytest.raises(table.TableError, match="panel markers"):
        table.load(shards["markers"], shards["manifest"])


def test_manifest_for_a_different_panel_fails(shards, tmp_path):
    other = tmp_path / "other.parquet"
    synth.make_panel(other, n_per_chrom=(180, 121))
    with pytest.raises(table.TableError, match="301"):
        table.load(other, shards["manifest"])


class FakeS3:
    """In-memory stand-in for table.S3Store: the same list/get contract, including
    FileNotFoundError for a key that is not there."""

    def __init__(self):
        self.objects = {}

    def put_dir(self, bucket, prefix, directory):
        for f in directory.iterdir():
            self.objects[(bucket, prefix + f.name)] = f.read_bytes()

    def list(self, bucket, prefix):
        return sorted(k for b, k in self.objects if b == bucket and k.startswith(prefix))

    def get(self, bucket, key):
        try:
            return self.objects[(bucket, key)]
        except KeyError:
            raise FileNotFoundError(f"s3://{bucket}/{key}") from None


@pytest.fixture
def s3(shards):
    store = FakeS3()
    store.put_dir("bkt", "results/SM/105943/markers/", shards["manifest"].parent)
    # manifests that are not this shard set's: a sibling prefix and a nested folder
    store.objects[("bkt", "results/SM/105943/markers-old/SM.markers.manifest.tsv")] = b"junk"
    store.objects[("bkt", "results/SM/105943/markers/archive/SM.markers.manifest.tsv")] = b"junk"
    return store


@pytest.mark.parametrize("uri", ["s3://bkt/results/SM/105943/markers/",
                                 "s3://bkt/results/SM/105943/markers",
                                 "s3://bkt/results/SM/105943/markers/SM.markers.manifest.tsv"])
def test_s3_shard_set_loads_the_same_sample_as_the_full_table(shards, s3, uri):
    d = table.load(shards["markers"], uri, s3=s3)
    assert llks(d) == llks(shards["full"])


def test_s3_missing_shard_fails_and_names_it(shards, s3):
    key = next(k for b, k in s3.objects if ".s004." in k)
    del s3.objects[("bkt", key)]
    with pytest.raises(table.TableError, match=key.rsplit("/", 1)[1]):
        table.load(shards["markers"], "s3://bkt/results/SM/105943/markers/", s3=s3)


def test_s3_prefix_without_a_manifest_fails(shards, s3):
    with pytest.raises(table.TableError, match="manifest"):
        table.load(shards["markers"], "s3://bkt/results/SM/999/markers/", s3=s3)


@pytest.mark.parametrize("n_manifests", [0, 2])
def test_directory_without_exactly_one_manifest_fails(shards, n_manifests):
    d = shards["manifest"].parent
    if n_manifests == 0:
        shards["manifest"].unlink()
    else:
        (d / "other.markers.manifest.tsv").write_text(shards["manifest"].read_text())
    with pytest.raises(table.TableError, match="manifest"):
        table.load(shards["markers"], d)
