"""remove_soft_clipped_bases_improved.py: soft clips and read-through past primer boundaries go."""
import remove_soft_clipped_bases_improved as scb

LEFTS, RIGHTS = [100, 300], [250, 450]


def sam(pos, cigar, seq, flag=0):
    return "\t".join(["r", str(flag), "chrM", str(pos), "60", cigar, "*", "0", "0", seq, "I" * len(seq)])


def run(pos, cigar, seq, flag=0, boundaries=(LEFTS, RIGHTS)):
    f = scb.process_sam_line(sam(pos, cigar, seq, flag), boundaries).split("\t")
    return int(f[3]), f[5], f[9]


def test_read_through_before_a_left_primer_is_cut():
    assert run(97, "20M", "AAA" + "C" * 17) == (100, "17M", "C" * 17)


def test_read_through_past_a_right_primer_is_cut():
    assert run(231, "22M", "C" * 20 + "GG") == (231, "20M", "C" * 20)


def test_reverse_reads_are_cut_by_position_too():
    assert run(98, "10M", "TT" + "C" * 8, flag=16) == (100, "8M", "C" * 8)


def test_reads_ending_on_boundaries_stay():
    assert run(100, "151M", "C" * 151) == (100, "151M", "C" * 151)


def test_a_product_spanning_two_amplicons_stays_whole():
    assert run(100, "351M", "C" * 351) == (100, "351M", "C" * 351)


def test_more_than_a_few_bases_outside_is_not_read_through():
    assert run(92, "20M", "C" * 20) == (92, "20M", "C" * 20)


def test_an_insertion_inside_the_cut_goes_with_it():
    assert run(98, "1M1I18M", "AGT" + "C" * 17) == (100, "17M", "C" * 17)


def test_soft_clips_go_first():
    assert run(98, "4S20M", "NNNN" + "AA" + "C" * 18) == (100, "18M", "C" * 18)


def test_without_primers_only_soft_clips_go():
    assert run(97, "4S20M", "NNNN" + "C" * 20, boundaries=None) == (97, "20M", "C" * 20)


def test_primer_boundaries_from_the_primer_files(reference_fasta):
    from pathlib import Path
    primers = Path(reference_fasta).parents[1] / "primers"
    lefts, rights = scb.primer_boundaries(reference_fasta, primers / "left_primers.fasta", primers / "right_primers_rc.fasta")
    assert 1 in lefts and 7819 in lefts and 7086 in lefts      # amplicon 001, 049, and 045's longer primer
    assert 182 in rights and 16622 in rights                   # amplicon 001, and 100 at the end of the extended reference
