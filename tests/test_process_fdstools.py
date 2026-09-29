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
mtNG_097 = chrM, 16094, 16276
[flanks]
"""


def run(tmp_path, rows, depth=10, tssv=None, reference_fasta=None, frame="separate"):
    """rows: (marker, sequence, reads) as in a sast.csv; tssv: (marker, sequence, reads)
    as in a tssv.csv, which switches the frames on. Returns the report table."""
    library = tmp_path / "library.txt"
    library.write_text(LIBRARY)
    sast = tmp_path / "sample.sast.csv"
    pd.DataFrame([{"marker": m, "sequence": s, "flags": "", "total": n, "total_mp_sum": 0}
                  for m, s, n in rows]).to_csv(sast, sep="\t", index=False)
    tssv_path = None
    if tssv is not None:
        tssv_path = tmp_path / "sample.tssv.csv"
        pd.DataFrame([{"marker": m, "sequence": s, "total": n, "forward": n, "reverse": 0}
                      for m, s, n in tssv]).to_csv(tssv_path, sep="\t", index=False)
    out = tmp_path / "sample.txt"
    fds.process_fdstools_sast(str(sast), str(library), str(out),
                              min_variant_frequency_pct=5.0, depth_threshold=depth,
                              length_heteroplasmy_threshold=10.0,
                              tssv_path=tssv_path and str(tssv_path),
                              reference_path=reference_fasta and str(reference_fasta), frame=frame)
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


# --- The frame regions (57-60, 300-315, 16180-16193) come from tssv.csv -------------

def amplicon(reference, start, end, edits):
    """rCRS start..end with edits {position: replacement}; '' deletes, 2+ bases insert."""
    return "".join(edits.get(p, reference[p - 1]) for p in range(start, end + 1))


def mtng_097(reference, edits):
    return amplicon(reference, 16094, 16276, edits)


@pytest.mark.parametrize("frame, expected", [
    ("separate", {"A16183-": 100.0, "T16189C": 100.0, "-16193.1C": 100.0, "T16217C": 100.0}),
    ("shared", {"A16183C": 100.0, "T16189C": 100.0, "T16217C": 100.0}),
])
def test_frame_rows_replace_fdstools_labels_in_the_region(tmp_path, reference, reference_fasta, frame, expected):
    # A3 C11 (spec example 4): the boundary shifted by one; T16217C lies outside the
    # region and keeps FDSTOOLS' label
    edits = {16183: "C", 16189: "C", 16217: "C"}
    report = run(tmp_path, [("mtNG_097", "A16183C T16189C T16217C", 500)],
                 tssv=[("mtNG_097", mtng_097(reference, edits), 500)], reference_fasta=reference_fasta, frame=frame)
    assert {k: v for k, v in calls(report).items() if k != "LOW"} == expected
    notes = report.set_index("FDSTOOLS")["variant_note"]
    assert notes["T16189C"] == f"{frame} frame" and pd.isna(notes["T16217C"])


def test_frame_rows_use_the_amplicon_reads_without_other_sequences(tmp_path, reference, reference_fasta):
    sast = [("mtNG_097", "T16189C", 70), ("mtNG_097", "T16189C 16193.1C", 30), ("mtNG_097", "Other sequences", 12)]
    tssv = [("mtNG_097", mtng_097(reference, {16189: "C"}), 70),
            ("mtNG_097", mtng_097(reference, {16189: "C", 16193: "CC"}), 30),
            ("mtNG_097", "Other sequences", 12)]
    report = run(tmp_path, sast, tssv=tssv, reference_fasta=reference_fasta)
    assert calls(report)["-16193.1c"] == 30.0
    assert report.set_index("FDSTOOLS").loc["-16193.1c", "rd_FDS"] == 30


def test_a_sample_with_changes_only_in_a_region_keeps_them(tmp_path, reference, reference_fasta):
    report = run(tmp_path, [("mtNG_003", "315.1C", 400)],
                 tssv=[("mtNG_003", amplicon(reference, 259, 367, {315: "CC"}), 400)], reference_fasta=reference_fasta)
    assert calls(report)["-315.1C"] == 100.0


def test_57_60_laid_from_the_left(tmp_path, reference, reference_fasta):
    # T57C plus one T more; the fewest-changes spelling would be -56.1C
    report = run(tmp_path, [("mtNG_001", "56.1C", 300)],
                 tssv=[("mtNG_001", amplicon(reference, 19, 155, {57: "CT"}), 300)], reference_fasta=reference_fasta)
    assert {k: v for k, v in calls(report).items() if k != "LOW"} == {"T57C": 100.0, "-60.1T": 100.0}
