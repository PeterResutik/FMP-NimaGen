"""score_mixtures.py: sites of a two-haplogroup mixture and how a label at a site is counted."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "simulation"))
import score_mixtures as sm  # noqa: E402

RCRS = sm.load_rcrs(sm.REPO / "resources/rCRS/rCRS_NimaGen.fasta")


def test_sites_are_single_base_differences_away_from_indels_and_c_stretches():
    major = sm.sequence("73G 263G 1000-", RCRS)
    minor = sm.sequence("263G 750G 1003G 310C", RCRS)
    assert [p for p, _, _ in sm.sites(major, minor, RCRS)] == [73, 750]


def test_the_minor_base_counts_as_found_only_as_the_iupac_code_of_both_bases():
    assert sm.status("A73R", "A", 73, "G", "A") == "minor found"
    assert sm.status("A73G", "A", 73, "G", "A") == "major only"
    assert sm.status(None, "A", 73, "A", "G") == "major only"
    assert sm.status("A73G", "A", 73, "A", "G") == "other"
