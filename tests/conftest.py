import os
import sys

import pytest

# Make microbiomeutils.py and the nrc package importable from the repository root.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Small feature table: 3 features x 2 samples, with a taxonomy column.
# The blank lines are deliberate (parser must ignore them).
FEATURE_TABLE = (
    "#FEATURE_ID\tS1\tS2\ttaxonomy\n"
    "\n"
    "f1\t1\t2\tk__Bacteria;p__Firmicutes;c__Bacilli\n"
    "f2\t3\t4\tk__Bacteria;p__Proteobacteria\n"
    "   \n"
    "f3\t0\t4\tk__Archaea;p__Euryarchaeota;c__Methanobacteria\n"
    "\n"
)

# Same counts, no taxonomy column.
FEATURE_TABLE_NO_TAX = (
    "#FEATURE_ID\tS1\tS2\n"
    "f1\t1\t2\n"
    "f2\t3\t4\n"
    "f3\t0\t4\n"
)


@pytest.fixture
def feature_table_file(tmp_path):
    p = tmp_path / "feature_table.tsv"
    p.write_text(FEATURE_TABLE)
    return p


@pytest.fixture
def feature_table_no_tax_file(tmp_path):
    p = tmp_path / "feature_table_notax.tsv"
    p.write_text(FEATURE_TABLE_NO_TAX)
    return p
