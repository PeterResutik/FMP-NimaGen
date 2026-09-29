"""process_fdstools_output_improved_better.py: FDSTOOLS sast rows -> report labels."""
import pandas as pd
import pytest

import process_fdstools_output_improved_better as fds

LIBRARY = """[genome_position]
mtNG_001 = chrM, 19, 155
mtNG_002 = chrM, 133, 267
mtNG_003 = chrM, 259, 367
mtNG_004 = chrM, 360, 480
mtNG_005 = chrM, 470, 610
mtNG_006 = chrM, 600, 800
[flanks]
"""


def run(tmp_path, rows, depth=10):
    """rows: (marker, sequence, reads) as in a sast.csv; returns the report table."""
    library = tmp_path / "library.txt"
    library.write_text(LIBRARY)
    sast = tmp_path / "sample.sast.csv"
    pd.DataFrame([{"marker": m, "sequence": s, "flags": "", "total": n, "total_mp_sum": 0}
                  for m, s, n in rows]).to_csv(sast, sep="\t", index=False)
    out = tmp_path / "sample.txt"
    fds.process_fdstools_sast(str(sast), str(library), str(out),
                              min_variant_frequency_pct=5.0, depth_threshold=depth,
                              length_heteroplasmy_threshold=10.0)
    return pd.read_csv(out, sep="\t")


def calls(report):
    return dict(zip(report["FDSTOOLS"], report["vf_FDS"]))


def low_amplicons(report):
    low = report[report["FDSTOOLS"] == "LOW"]
    return dict(zip(low["marker"], low["rd_FDS"]))


@pytest.mark.parametrize("threshold", [10, 90])
def test_lh_bounds_same_for_either_side(threshold):
    assert fds.lh_bounds_pct(threshold) == (10, 90)


def test_frequency_counts_reads_without_other_sequences(tmp_path):
    report = run(tmp_path, [("mtNG_001", "T152C", 498), ("mtNG_001", "Other sequences", 27)])
    assert calls(report)["T152C"] == 100.0


def test_low_counts_all_rows_of_an_amplicon_and_keeps_its_calls(tmp_path):
    # 6 reads plus an "Other sequences" row: flagged, and the call is still reported
    report = run(tmp_path, [("mtNG_002", "T204C", 6), ("mtNG_002", "Other sequences", 20)])
    assert low_amplicons(report)["mtNG_002"] == 6
    assert calls(report)["T204C"] == 100.0


def test_amplicon_missing_from_the_file_is_low(tmp_path):
    report = run(tmp_path, [("mtNG_001", "T152C", 498)])
    assert low_amplicons(report)["mtNG_004"] == 0
    assert "mtNG_001" not in low_amplicons(report)


@pytest.mark.parametrize("with_insertion, expected", [
    (5, {}),
    (50, {"-315.1c": 50.0}),
    (95, {"-315.1C": 95.0}),
])
def test_insertion_floor_minor_major(tmp_path, with_insertion, expected):
    report = run(tmp_path, [("mtNG_003", "A263G 315.1C", with_insertion),
                            ("mtNG_003", "A263G", 100 - with_insertion)])
    length_calls = {k: v for k, v in calls(report).items() if "." in k}
    assert length_calls == expected


@pytest.mark.parametrize("with_deletion, expected", [
    (5, {}),
    (40, {"A523a": 40.0, "C524c": 40.0}),
    (100, {"A523-": 100.0, "C524-": 100.0}),
])
def test_deletion_floor_minor_major(tmp_path, with_deletion, expected):
    rows = [("mtNG_005", "A523DEL C524DEL", with_deletion)]
    if with_deletion < 100:
        rows.append(("mtNG_005", "REF", 100 - with_deletion))
    report = run(tmp_path, rows)
    assert {k: v for k, v in calls(report).items() if k != "LOW"} == expected


@pytest.mark.parametrize("with_alt, expected", [(29, "C756M"), (97, "C756A")])
def test_substitution_minor_iupac_major_base(tmp_path, with_alt, expected):
    report = run(tmp_path, [("mtNG_006", "C756A", with_alt), ("mtNG_006", "REF", 100 - with_alt)])
    assert calls(report)[expected] == with_alt
