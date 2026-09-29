"""notation.py: how one molecule differs from rCRS (the general rule)."""
import json
from pathlib import Path

import pytest

import notation

WINDOWS = json.loads((Path(__file__).parent / "fixtures" / "notation_mitoleaf_windows.json").read_text())


@pytest.fixture
def rcrs(reference):
    return "".join(reference)


def molecule(rcrs, start, end, edits):
    """rCRS start..end with edits {position: replacement}; '' deletes a base,
    two or more bases insert after the first. The N at 3107 has no base unless edited."""
    out = ""
    for p in range(start, end + 1):
        base = rcrs[p - 1]
        out += edits.get(p, "" if base == "N" else base)
    return out


@pytest.mark.parametrize("start, end, edits, expected", [
    (60, 85, {73: "G"}, ["A73G"]),
    (560, 585, {573: "CC"}, ["-573.1C"]),                                   # insertion at the 3' end of the C-run
    (505, 535, {515: "", 516: ""}, ["A523-", "C524-"]),                     # AC unit off the 3'-most copy
    (505, 535, {524: "CAC"}, ["-524.1A", "-524.2C"]),
    (8265, 8300, {p: "" for p in range(8272, 8281)},                        # 9-bp repeat copy, 3'-most copy
     ["C8281-", "C8282-", "C8283-", "C8284-", "C8285-", "T8286-", "C8287-", "T8288-", "A8289-"]),
    (15925, 15955, {15940: "C", 15943: ""}, ["T15940C", "T15944-"]),        # L3f3: fewest changes
    (70, 90, {80: "G", 81: "C"}, ["C80G", "G81C"]),                         # tie: substitutions beat gaps
])
def test_general_rule(rcrs, start, end, edits, expected):
    assert notation.describe(molecule(rcrs, start, end, edits), rcrs, start, end) == expected


@pytest.mark.parametrize("base, expected", [("T", ["-3109.1T"]), ("C", ["-3106.1C"])])
def test_a_base_at_the_rcrs_n_is_an_insertion(rcrs, base, expected):
    # L3y: its extra T sits where rCRS has the N placeholder
    assert notation.describe(molecule(rcrs, 3095, 3120, {3107: base}), rcrs, 3095, 3120) == expected


@pytest.mark.parametrize("edits, expected", [
    ({452: ""}, ["T455-"]),                 # one T fewer
    ({452: "TT"}, ["-455.1T"]),             # one T more
    ({455: "C"}, ["T455C"]),                # boundary shift: a substitution in the shared frame
    ({460: "C"}, ["T460C"]),
    ({452: "", 460: "C"}, ["T455-", "T460C"]),
])
def test_452_463_leading_t_run(rcrs, edits, expected):
    assert notation.describe(molecule(rcrs, 445, 470, edits), rcrs, 445, 470) == expected


def test_nothing_is_shifted_across_the_origin(rcrs):
    # 16569 and 1 are both G; a deleted G stays at 16569
    assert notation.describe(molecule(rcrs, 16550, 16569, {16569: ""}), rcrs, 16550, 16569) == ["G16569-"]


def test_rcrs_has_no_changes(rcrs):
    assert notation.describe(molecule(rcrs, 8265, 8300, {}), rcrs, 8265, 8300) == []


@pytest.mark.parametrize("case", WINDOWS, ids=lambda c: f"{c['example']}:{c['start']}-{c['end']}")
def test_mitoleaf_windows(rcrs, case):
    """Every distinct mitoLEAF window with a length change outside 16180-16193,
    300-315 and 57-60. The expected labels equal mitoLEAF's except in 10 windows,
    where mitoLEAF follows the order of changes in the tree (G247A A249- for G247-)."""
    assert notation.describe(case["molecule"], rcrs, case["start"], case["end"]) == case["expected"]
