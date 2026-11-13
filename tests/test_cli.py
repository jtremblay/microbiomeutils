"""Tests for the microbiomeutils.py command-line subcommands."""
import io

import numpy as np
import pytest
from skbio import DistanceMatrix, OrdinationResults
from skbio.io import LSMatFormatError

import microbiomeutils as mu


def run_main(argv, capsys):
    mu.main([str(a) for a in argv])
    return capsys.readouterr().out


# ---------------------------------------------------------------- argument errors

@pytest.mark.parametrize("argv", [
    ["betadiv"],
    ["taxsum"],
    ["pcoa"],
])
def test_missing_input_is_usage_error(argv):
    with pytest.raises(SystemExit) as exc:
        mu.main(argv)
    assert exc.value.code == 2


@pytest.mark.parametrize("metric", ["weighted-unifrac", "unweighted-unifrac"])
def test_unifrac_without_tree_is_usage_error(metric, feature_table_file, capsys):
    with pytest.raises(SystemExit) as exc:
        mu.main(["betadiv", "-i", str(feature_table_file), "-m", metric])
    assert exc.value.code == 2
    assert "--infile-tree" in capsys.readouterr().err


def test_emperor_requires_outdir(feature_table_file):
    with pytest.raises(SystemExit) as exc:
        mu.main(["emperor", "-i", str(feature_table_file), "-m", str(feature_table_file)])
    assert exc.value.code == 2


# ---------------------------------------------------------------- taxsum

def test_taxsum_absolute(feature_table_file, capsys):
    out = run_main(["taxsum", "-i", feature_table_file, "-l", "1"], capsys)
    lines = [line.split("\t") for line in out.splitlines()]
    assert lines[0] == ["Taxon", "S1", "S2"]
    rows = {r[0]: [float(x) for x in r[1:]] for r in lines[1:]}
    assert rows == {"k__Archaea": [0.0, 4.0], "k__Bacteria": [4.0, 6.0]}


def test_taxsum_relative(feature_table_file, capsys):
    out = run_main(["taxsum", "-i", feature_table_file, "-l", "1", "-t", "relative"], capsys)
    rows = {r[0]: [float(x) for x in r[1:]] for r in (l.split("\t") for l in out.splitlines()[1:])}
    assert rows["k__Bacteria"] == [1.0, 0.6]
    assert rows["k__Archaea"] == [0.0, 0.4]


def test_taxsum_without_taxonomy_raises(feature_table_no_tax_file):
    with pytest.raises(ValueError, match="No taxonomy"):
        mu.main(["taxsum", "-i", str(feature_table_no_tax_file)])


# ---------------------------------------------------------------- betadiv / pcoa

def test_betadiv_bray_curtis(feature_table_file, capsys):
    out = run_main(["betadiv", "-i", feature_table_file], capsys)
    dm = DistanceMatrix.read(io.StringIO(out))
    assert dm.ids == ("S1", "S2")
    # S1 = [1, 3, 0], S2 = [2, 4, 4]: sum|diff| / sum = (1 + 1 + 4) / 14
    assert dm["S1", "S2"] == pytest.approx(6 / 14)


def test_betadiv_unweighted_unifrac(tmp_path, capsys):
    table = tmp_path / "t.tsv"
    table.write_text("#FEATURE_ID\tS1\tS2\nf1\t5\t0\nf2\t0\t5\nf3\t0\t0\n")
    tree = tmp_path / "tree.nwk"
    tree.write_text("((f1:1,f2:1):1,f3:2);\n")
    out = run_main(["betadiv", "-i", table, "-m", "unweighted-unifrac", "-t", tree], capsys)
    dm = DistanceMatrix.read(io.StringIO(out))
    # Observed branches: f1 (1), f2 (1), shared internal (1); unique = 2 of 3.
    assert dm["S1", "S2"] == pytest.approx(2 / 3)


def test_pcoa_roundtrip(tmp_path, capsys):
    ids = ["A", "B", "C", "D"]
    data = np.array([
        [0.0, 0.2, 0.6, 0.7],
        [0.2, 0.0, 0.5, 0.6],
        [0.6, 0.5, 0.0, 0.3],
        [0.7, 0.6, 0.3, 0.0],
    ])
    dm_file = tmp_path / "dm.tsv"
    DistanceMatrix(data, ids=ids).write(str(dm_file))
    out = run_main(["pcoa", "-i", dm_file], capsys)
    ord_res = OrdinationResults.read(io.StringIO(out))
    assert list(ord_res.samples.index) == ids
    assert ord_res.proportion_explained.iloc[0] > 0.5


def test_pcoa_rejects_mismatched_row_labels(tmp_path):
    dm_file = tmp_path / "bad.tsv"
    dm_file.write_text("\tA\tB\nB\t0.0\t0.5\nA\t0.5\t0.0\n")
    with pytest.raises(LSMatFormatError, match="mismatched IDs"):
        mu.main(["pcoa", "-i", str(dm_file)])
