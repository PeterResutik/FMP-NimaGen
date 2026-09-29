"""frames.py: placing a molecule (an amplicon's sequence from FDSTOOLS) into a region."""
import json
from pathlib import Path

import pytest

import frames
import merge_fdstools_mutect2_improved as merge

HERE = Path(__file__).parent
AMPLICONS = json.loads((HERE / "fixtures" / "placement_mitoleaf_amplicons.json").read_text())
LIBRARY = HERE.parent / "resources" / "fdstools" / "mtNG_library_file.txt"
R16189, R310, R57 = frames.REGIONS["16180-16193"], frames.REGIONS["300-315"], frames.REGIONS["57-60"]


@pytest.fixture
def rcrs(reference):
    return "".join(reference)


def amplicon(rcrs, start, end, edits={}):
    """rCRS start..end with edits {position: replacement}; '' deletes, 2+ bases insert after the first."""
    return "".join(edits.get(p, rcrs[p - 1]) for p in range(start, end + 1))


def test_each_region_lies_in_one_amplicon():
    ranges = merge.load_amplicon_ranges(LIBRARY)
    assert frames.covering(R57, ranges) == ["mtNG_001"]
    assert frames.covering(R310, ranges) == ["mtNG_003"]
    assert frames.covering(R16189, ranges) == ["mtNG_097"]


@pytest.mark.parametrize("edits, expected", [
    ({}, "AAAACCCCCTCCCC"),
    ({16177: "G"}, "AAAACCCCCTCCCC"),           # variant in the flank CACATC the old cut relied on
    ({16197: "A"}, "AAAACCCCCTCCCC"),           # ... and in ATGCTT on the other side
    ({16179: "CA"}, "AAAAACCCCCTCCCC"),         # an A inserted before the region joins its A-run
    ({16179: "CC"}, "AAAACCCCCTCCCC"),          # a C there lengthens the run in front, not the region
    ({16179: ""}, "AAAACCCCCTCCCC"),            # deletion right before the region
    ({16193: "CC"}, "AAAACCCCCTCCCCC"),         # insertion right after the region belongs to it
    ({16120: "", 16121: ""}, "AAAACCCCCTCCCC"),  # a deletion upstream moves every later base
    ({16110: "TAAA"}, "AAAACCCCCTCCCC"),        # ... and so does an insertion
    ({16183: "C", 16189: "C", 16193: "CC", 16217: "C"}, "AAACCCCCCCCCCCC"),
])
def test_16180_16193(rcrs, edits, expected):
    assert frames.region_bases(amplicon(rcrs, 16094, 16276, edits), 16094, 16276, R16189, rcrs) == expected


def test_300_315_and_57_60(rcrs):
    assert frames.region_bases(amplicon(rcrs, 259, 367, {315: "CC"}), 259, 367, R310, rcrs) == "AAACCCCCCCTCCCCCC"
    assert frames.region_bases(amplicon(rcrs, 19, 155, {57: "CT"}), 19, 155, R57, rcrs) == "CTTTT"


def test_molecules_are_summed_per_region_sequence(rcrs):
    sequences = [
        (amplicon(rcrs, 16094, 16276, {16189: "C"}), 60),
        (amplicon(rcrs, 16094, 16276, {16189: "C", 16217: "C"}), 30),   # differs outside the region only
        (amplicon(rcrs, 16094, 16276), 10),
        ("Other sequences", 7),
    ]
    assert frames.region_molecules(sequences, 16094, 16276, R16189, rcrs) == [
        ("AAAACCCCCCCCCC", 90), ("AAAACCCCCTCCCC", 10)]


@pytest.mark.parametrize("case", AMPLICONS, ids=lambda c: f"{c['region']}:{c['example']}")
def test_mitoleaf_amplicons(rcrs, case):
    """The whole amplicon of one example haplogroup for every distinct mitoLEAF
    molecule in the three regions; placement must give that molecule."""
    region = frames.REGIONS[case["region"]]
    assert frames.region_bases(case["sequence"], case["start"], case["end"], region, rcrs) == case["molecule"]
