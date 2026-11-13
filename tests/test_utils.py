"""Unit tests for the parsing and taxonomic-summary helpers in nrc/utils.py."""
import io

import numpy as np
import pytest

from conftest import FEATURE_TABLE, FEATURE_TABLE_NO_TAX
from nrc.utils import (
    FeatureMissingError,
    convert_feature_table_relative,
    getFeatureTable,
    make_summary,
    parse_feature_table,
    sum_counts_by_consensus,
)


def test_parse_feature_table_with_taxonomy_skips_blank_lines():
    sample_ids, feature_ids, counts, lineages = parse_feature_table(io.StringIO(FEATURE_TABLE))
    assert sample_ids == ["S1", "S2"]
    assert feature_ids == ["f1", "f2", "f3"]
    assert np.array(counts).tolist() == [[1, 2], [3, 4], [0, 4]]
    assert lineages[1] == ["k__Bacteria", "p__Proteobacteria"]


def test_parse_feature_table_without_taxonomy():
    sample_ids, feature_ids, counts, lineages = parse_feature_table(io.StringIO(FEATURE_TABLE_NO_TAX))
    assert sample_ids == ["S1", "S2"]
    assert feature_ids == ["f1", "f2", "f3"]
    assert np.array(counts).tolist() == [[1, 2], [3, 4], [0, 4]]
    assert lineages == []


@pytest.mark.parametrize("header", ["#FEATURE_ID", "#OTU ID", "#FEATURE ID", "#FEATUREID"])
def test_parse_feature_table_accepts_header_variants(header):
    table = FEATURE_TABLE_NO_TAX.replace("#FEATURE_ID", header, 1)
    sample_ids, feature_ids, _, _ = parse_feature_table(io.StringIO(table))
    assert sample_ids == ["S1", "S2"]
    assert feature_ids == ["f1", "f2", "f3"]


def test_parse_feature_table_float_counts():
    table = "#FEATURE_ID\tS1\tS2\nf1\t0.5\t1.5\n"
    _, _, counts, _ = parse_feature_table(io.StringIO(table))
    assert np.array(counts).tolist() == [[0.5, 1.5]]


def test_get_feature_table_missing_file_raises():
    with pytest.raises(FeatureMissingError):
        getFeatureTable("/nonexistent/path/feature_table.tsv")


def test_sum_counts_by_consensus_truncates_and_pads():
    table = parse_feature_table(io.StringIO(FEATURE_TABLE))
    # Lineages can be summarised more than once (they used to be one-shot iterators).
    sum_counts_by_consensus(table, 1)
    result, _ = sum_counts_by_consensus(table, 3)
    assert result[("k__Bacteria", "p__Firmicutes", "c__Bacilli")].tolist() == [1, 2]
    # Lineage shorter than the requested level is padded with 'Other'.
    assert result[("k__Bacteria", "p__Proteobacteria", "Other")].tolist() == [3, 4]


def test_make_summary_level1():
    table = parse_feature_table(io.StringIO(FEATURE_TABLE))
    summary, header = make_summary(table, 1, None, None)
    assert header == ["Taxon", "S1", "S2"]
    rows = {row[0]: list(row[1:]) for row in summary}
    assert rows == {("k__Archaea",): [0, 4], ("k__Bacteria",): [4, 6]}


def test_convert_feature_table_relative_columns_sum_to_one():
    table = parse_feature_table(io.StringIO(FEATURE_TABLE))
    _, _, rel, _ = convert_feature_table_relative(table)
    np.testing.assert_allclose(rel.sum(axis=0), [1.0, 1.0])
    np.testing.assert_allclose(rel[:, 0], [0.25, 0.75, 0.0])


def test_convert_feature_table_relative_empty_sample_is_zero():
    table = "#FEATURE_ID\tS1\tS2\nf1\t0\t2\nf2\t0\t2\n"
    _, _, rel, _ = convert_feature_table_relative(parse_feature_table(io.StringIO(table)))
    assert rel[:, 0].tolist() == [0.0, 0.0]  # no NaN from 0/0
    np.testing.assert_allclose(rel[:, 1], [0.5, 0.5])
