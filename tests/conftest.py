"""Shared fixtures. pytest.ini puts resources/scripts on the import path."""
from pathlib import Path

import pytest
from Bio import SeqIO

REFERENCE_FASTA = Path(__file__).resolve().parents[1] / "resources" / "rCRS" / "rCRS_NimaGen.fasta"


@pytest.fixture
def reference_fasta():
    return REFERENCE_FASTA


@pytest.fixture
def reference():
    """rCRS_NimaGen as a list of bases (index = position - 1). A fresh copy per
    test, because the report functions write substitutions into it."""
    return list(str(SeqIO.read(REFERENCE_FASTA, "fasta").seq))
