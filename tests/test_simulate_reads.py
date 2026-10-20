"""simulate_reads.py: reads with the qualities and errors of a real read (--error-profile)."""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "simulation"))
import simulate_reads as sr  # noqa: E402


def test_the_read_takes_the_qualities_and_errors_of_the_real_read():
    read, qual = sr.sequence_read_like("ACGTACGT", random.Random(1), ("IIII#III", [(1, "T"), (4, "G")]))
    assert read == "ATGTGCGT"
    assert qual == "IIII#III"


def test_an_error_towards_the_base_the_molecule_has_becomes_another_base():
    for seed in range(20):
        read, _ = sr.sequence_read_like("ACGT", random.Random(seed), ("IIII", [(2, "G")]))
        assert read[2] not in "G" and read[:2] + read[3:] == "ACT"


def test_a_molecule_longer_than_the_real_read_repeats_its_last_quality():
    read, qual = sr.sequence_read_like("ACGTAC", random.Random(1), ("I#", [(5, "A")]))
    assert qual == "I#####"
    assert read == "ACGTAA"


def test_errors_beyond_the_molecule_are_ignored():
    read, _ = sr.sequence_read_like("ACG", random.Random(1), ("IIIIII", [(4, "T")]))
    assert read == "ACG"


def test_a_deletion_skips_molecule_bases_and_the_read_reads_on():
    read, qual = sr.sequence_read_like("ACGTACGTAA", random.Random(1), ("IIIIII", [], [(2, "D", 1)]), 6)
    assert read == "ACTACG"
    assert qual == "IIIIII"


def test_an_insertion_takes_the_real_reads_inserted_bases():
    read, _ = sr.sequence_read_like("ACGTACGT", random.Random(1), ("IIIIII", [], [(2, "I", "TT")]), 6)
    assert read == "ACTTGT"


def test_a_read_with_a_deletion_still_ends_at_the_molecule_end():
    read, qual = sr.sequence_read_like("ACGT", random.Random(1), ("I#II", [], [(1, "D", 2)]), 4)
    assert read == "AT"
    assert qual == "I#"


def test_substitution_errors_stay_at_their_cycle_after_an_insertion():
    read, _ = sr.sequence_read_like("CCCCCC", random.Random(1), ("IIIIII", [(4, "A")], [(1, "I", "G")]), 6)
    assert read == "CGCCAC"


def test_a_profile_without_insertions_and_deletions_reads_up_to_the_length():
    read, qual = sr.sequence_read_like("ACGTAC", random.Random(1), ("II", [(1, "T")]), 4)
    assert read == "ATGT"
    assert qual == "IIII"
