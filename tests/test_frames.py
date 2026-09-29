"""frames.py: the shared and the separate frame in 16180-16193, 300-315 and 57-60."""
import json
from pathlib import Path

import pytest

import frames

MOLECULES = json.loads((Path(__file__).parent / "fixtures" / "frame_mitoleaf_molecules.json").read_text())
R16189, R310, R57 = frames.REGIONS["16180-16193"], frames.REGIONS["300-315"], frames.REGIONS["57-60"]


@pytest.fixture
def rcrs(reference):
    return "".join(reference)


def mol(a, c, x=None, c2=None):
    """A leading A-run, a C-run and optionally an interrupt base with a second C-run."""
    return "A" * a + "C" * c + ("" if x is None else x + "C" * c2)


def report(region, molecules, frame, rcrs):
    """Rows as the spec writes them: majors plain, minors with their rounded share."""
    coverage = sum(reads for _, reads in molecules)
    return [label if label[-1] in "ACGT-" and share >= 90 else f"{label} {share:.0f}%"
            for label, share in frames.rows(region, molecules, coverage, frame, rcrs)]


# The worked examples of the spec (16180-16193); reads out of 100
SPEC = [
    (1, [(mol(4, 5, "T", 4), 100)], [], []),
    (2, [(mol(4, 10), 100)], ["T16189C"], ["T16189C"]),
    (3, [(mol(4, 11), 100)], ["T16189C", "-16193.1C"], ["T16189C", "-16193.1C"]),
    (4, [(mol(3, 11), 100)], ["A16183C", "T16189C"], ["A16183-", "T16189C", "-16193.1C"]),
    (5, [(mol(4, 10), 70), (mol(4, 11), 30)], ["T16189C", "-16193.1c 30%"], ["T16189C", "-16193.1c 30%"]),
    (6, [(mol(3, 11), 70), (mol(4, 10), 30)], ["A16183M 70%", "T16189C"], ["A16183a 70%", "T16189C", "-16193.1c 70%"]),
    (7, [(mol(4, 5, "T", 4), 60), (mol(4, 10), 40)], ["T16189Y 40%"], ["T16189Y 40%"]),
    (8, [(mol(2, 10), 10), (mol(3, 10), 20), (mol(3, 11), 30), (mol(2, 12), 20), (mol(3, 12), 10), (mol(3, 13), 10)],
     ["A16182M 30%", "A16183C", "T16189C", "C16192c 10%", "C16193c 30%", "-16193.1c 20%", "-16193.2c 10%"],
     ["A16182a 30%", "A16183-", "T16189C", "-16193.1c 70%", "-16193.2c 40%", "-16193.3c 10%"]),
    (9, [("AAACA" + "CCCCCTCCCC", 100)], ["A16183C", "C16184A", "-16188.1C"], ["A16183C", "C16184A", "-16188.1C"]),
    (10, [(mol(3, 6, "T", 4), 100)], ["A16183C"], ["A16183-", "-16188.1C"]),
    (11, [(mol(3, 5, "T", 4), 100)], ["A16183C", "C16188-"], ["A16183-"]),
    (12, [(mol(4, 5, "T", 5), 100)], ["-16193.1C"], ["-16193.1C"]),
    (13, [(mol(4, 5, "T", 4), 60), (mol(4, 10), 25), (mol(4, 11), 10), (mol(4, 12), 5)],
     ["T16189Y 40%", "-16193.1c 15%"], ["T16189Y 40%", "-16193.1c 15%"]),
    # three bases at 16189: one row per base until the three-base IUPAC code (spec point 9)
    (14, [(mol(4, 5, "T", 4), 40), (mol(4, 10), 30), (mol(3, 6, "A", 4), 10), (mol(4, 11), 10), (mol(4, 5, "T", 5), 10)],
     ["A16183M 10%", "T16189Y 40%", "T16189W 10%", "-16193.1c 20%"],
     ["A16183a 10%", "-16188.1c 10%", "T16189Y 40%", "T16189W 10%", "-16193.1c 20%"]),
    (15, [("AAACACCCCCCCCC", 63), ("AAACACCCCCCCCCC", 30), ("AAACACCCCCCCCCCC", 4), ("AAACACCCCCCCC", 4)],
     ["A16183C", "C16184A", "T16189C", "-16193.1c 34%"], ["A16183C", "C16184A", "T16189C", "-16193.1c 34%"]),
    (16, [("AAAACCCTCCCCCC", 100)], ["C16187T", "T16189C"], ["C16187T", "T16189C"]),
]


@pytest.mark.parametrize("example, molecules, shared, separate", SPEC, ids=[f"example{s[0]}" for s in SPEC])
def test_spec_examples(rcrs, example, molecules, shared, separate):
    assert report(R16189, molecules, "shared", rcrs) == shared
    assert report(R16189, molecules, "separate", rcrs) == separate


@pytest.mark.parametrize("molecule, shared, separate", [
    (mol(3, 7, "T", 6), ["-315.1C"], ["-315.1C"]),
    (mol(3, 8, "T", 6), ["-309.1C", "-315.1C"], ["-309.1C", "-315.1C"]),
    (mol(3, 6, "T", 5), ["C309-"], ["C309-"]),
    (mol(3, 13), ["T310C"], ["T310C"]),
    (mol(2, 8, "T", 5), ["A302C"], ["A302-", "-309.1C"]),       # boundary shift
])
def test_300_315(rcrs, molecule, shared, separate):
    assert report(R310, [(molecule, 100)], "shared", rcrs) == shared
    assert report(R310, [(molecule, 100)], "separate", rcrs) == separate


@pytest.mark.parametrize("shift, shared, separate", [
    (7, ["A302M 7%", "-315.1C"], ["-315.1C"]),
    (12, ["A302M 12%", "-315.1C"], ["A302a 12%", "-309.1c 12%", "-315.1C"]),
])
def test_boundary_shift_uses_each_spellings_threshold(rcrs, shift, shared, separate):
    # a substitution in the shared frame (from --min_vf, 5%), a length change in the
    # separate frame (from --lh_thresh, 10%); spec point 15 as revised
    molecules = [(mol(3, 7, "T", 6), 100 - shift), (mol(2, 8, "T", 6), shift)]
    assert report(R310, molecules, "shared", rcrs) == shared
    assert report(R310, molecules, "separate", rcrs) == separate


@pytest.mark.parametrize("molecule, expected", [("CTTTT", ["T57C", "-60.1T"]), ("TTT", ["T60-"])])
def test_57_60_laid_from_the_left_in_both_frames(rcrs, molecule, expected):
    assert report(R57, [(molecule, 100)], "shared", rcrs) == expected
    assert report(R57, [(molecule, 100)], "separate", rcrs) == expected


@pytest.mark.parametrize("molecule, expected", [
    ("AAAACCCCCTCCTC", ["C16192T"]),              # T kept, a second non-C base
    ("AAAACCCCCCCCTC", ["T16189C", "C16192T"]),   # L2a1: laid from the left, not anchored at the T
])
def test_anchor_only_where_it_saves_changes(rcrs, molecule, expected):
    for frame in ("shared", "separate"):
        assert report(R16189, [(molecule, 100)], frame, rcrs) == expected


def test_on_a_tie_only_the_separate_frame_anchors(rcrs):
    # M13a: one C fewer before the T, one more after it; both readings need two changes
    molecule = "AAAACCCCTCCCCC"
    assert report(R16189, [(molecule, 100)], "shared", rcrs) == ["C16188T", "T16189C"]    # as mitoLEAF
    assert report(R16189, [(molecule, 100)], "separate", rcrs) == ["C16188-", "-16193.1C"]  # the T stays


def test_of_two_equally_close_candidates_the_cheaper_anchors(rcrs):
    # Z1a1a: in the separate reading both T's are two places from 16189
    assert report(R16189, [("AACCCTCCCTCCCC", 100)], "separate", rcrs) == [
        "A16182-", "A16183-", "C16187T", "-16188.1C", "-16188.2C"]


def test_s26_02989(rcrs):
    """Eight molecules without the T, the largest at 31% (a real sample)."""
    molecules = [(mol(2, 12), 33), (mol(2, 11), 22), (mol(3, 11), 15), (mol(2, 13), 9),
                 (mol(3, 12), 9), (mol(2, 10), 8), (mol(3, 10), 7), (mol(3, 13), 3)]
    assert report(R16189, molecules, "shared", rcrs) == [
        "A16182M 68%", "A16183C", "T16189C", "C16193c 35%", "-16193.1c 20%"]
    assert report(R16189, molecules, "separate", rcrs) == [
        "A16182a 68%", "A16183-", "T16189C", "-16193.1c 86%", "-16193.2c 51%", "-16193.3c 11%"]


@pytest.mark.parametrize("case", MOLECULES, ids=lambda c: f"{c['region']}:{c['example']}")
def test_mitoleaf_molecules(rcrs, case):
    """Every distinct molecule the mitoLEAF haplogroups have in these regions. The
    shared frame equals mitoLEAF's labels except for L5a (mitoLEAF writes A16183-),
    U2e3a (a tie read the other way) and M39 in 57-60."""
    region = frames.REGIONS[case["region"]]
    assert frames.labels(region, case["molecule"], "shared", rcrs) == case["shared"]
    assert frames.labels(region, case["molecule"], "separate", rcrs) == case["separate"]
