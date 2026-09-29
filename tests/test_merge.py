"""merge_fdstools_mutect2_improved.py: both callers' tables -> one report."""
import subprocess
import sys

import pandas as pd
import pytest
from openpyxl import load_workbook

import merge_fdstools_mutect2_improved as merge

AMPLICONS = {"mtNG_001": (19, 155), "mtNG_002": (133, 267), "mtNG_003": (259, 367),
             "mtNG_005": (470, 610), "mtNG_097": (16094, 16276)}
MT2_COLUMNS = ["MUTECT2", "vf_MT2", "rd_MT2", "MBQ", "Filter", "Pos", "Ref", "Variant", "GT", "Type"]


def fds(label, vf, depth, marker):
    """A row of the FDSTOOLS report table (vf in percent)."""
    start, end = AMPLICONS[marker]
    return {"FDSTOOLS": label, "vf_FDS": vf, "rd_FDS": round(vf * depth / 100),
            "interpolated_total_coverage": depth, "marker": marker,
            "marker_range": f"{start}-{end}", "num_markers": 1}


def fds_low(marker, reads):
    start, end = AMPLICONS[marker]
    return {"FDSTOOLS": "LOW", "rd_FDS": reads, "marker": marker, "marker_range": f"{start}-{end}"}


def mt2(label, vf, ad, pos, ref, alt, kind):
    """A row of the Mutect2 report table (vf as a fraction, ad as ref,alt reads)."""
    return dict(zip(MT2_COLUMNS, [label, vf, ad, "35,35", "PASS", pos, ref, alt, "0/1", kind]))


def write_inputs(tmp_path, fds_rows, mt2_rows, depths):
    library = tmp_path / "library.txt"
    library.write_text("[genome_position]\n" + "".join(
        f"{m} = chrM, {s}, {e}\n" for m, (s, e) in AMPLICONS.items()) + "[flanks]\n")
    fds_file, mt2_file, depth_file = tmp_path / "fds.txt", tmp_path / "mt2.txt", tmp_path / "depth.txt"
    pd.DataFrame(fds_rows).to_csv(fds_file, sep="\t", index=False)
    pd.DataFrame(mt2_rows, columns=MT2_COLUMNS).to_csv(mt2_file, sep="\t", index=False)
    # p09's depth at each amplicon's middle position; 500 reads unless given
    depth_file.write_text("".join(f"chrM\t{(s + e) // 2}\t{depths.get(m, 500)}\n"
                                  for m, (s, e) in AMPLICONS.items()))
    return fds_file, mt2_file, depth_file, library


def run(tmp_path, fds_rows, mt2_rows, depths={}, weighted=False):
    fds_file, mt2_file, depth_file, library = write_inputs(tmp_path, fds_rows, mt2_rows, depths)
    return merge.merge_variant_callers(str(fds_file), str(mt2_file), lh_thresh=10.0, min_vf=5.0,
                                       mutect2_depth_file=str(depth_file), marker_map=str(library),
                                       depth_threshold=10, weighted_average=weighted)


def calls(report):
    return report[report["FMP"] != "LOW"].set_index("FMP")


@pytest.mark.parametrize("a, b", [("-309.1C", "-309.1c"), ("T16189C", "T16189Y"), ("A523-", "A523a")])
def test_merge_key_same_locus_any_spelling(a, b):
    assert merge.merge_key(a) == merge.merge_key(b)


def test_merge_key_keeps_different_alleles_apart():
    assert merge.merge_key("C756M") != merge.merge_key("C756Y")


def test_agreeing_callers_share_a_row(tmp_path):
    report = calls(run(tmp_path, [fds("T152C", 100.0, 498, "mtNG_001")],
                       [mt2("T152C", 0.998, "1,497", 152, "T", "C", "SNP")]))
    assert len(report) == 1
    assert report.loc["T152C", "called_by_FDSTOOLS"] is True
    assert report.loc["T152C", "called_by_MUTECT2"] is True


def test_silent_caller_is_false(tmp_path):
    report = calls(run(tmp_path, [fds("A263G", 100.0, 900, "mtNG_003")], []))
    assert report.loc["A263G", "called_by_MUTECT2"] is False


def test_length_disagreement_decided_by_average(tmp_path):
    # FDSTOOLS 92% (major), Mutect2 80% (minor): average 86% < 90% -> minor
    report = calls(run(tmp_path, [fds("-309.1C", 92.0, 900, "mtNG_003")],
                       [mt2("-309.1c", 0.80, "180,720", 302, "A", "AC", "LHP")]))
    assert list(report.index) == ["-309.1c"]
    assert report.loc["-309.1c", "Type"] == "LHP"
    assert report.loc["-309.1c", "called_by_FDSTOOLS"] == "DISAGREEMENT"
    assert report.loc["-309.1c", "called_by_MUTECT2"] is True


@pytest.mark.parametrize("weighted, reported, kind, disagreeing", [
    (False, "T16189Y", "PHP", "called_by_FDSTOOLS"),   # (100 + 66.7) / 2 < 95
    (True, "T16189C", "SNP", "called_by_MUTECT2"),     # 22 FDSTOOLS reads outweigh 1 Mutect2 read
])
def test_substitution_disagreement(tmp_path, weighted, reported, kind, disagreeing):
    report = calls(run(tmp_path, [fds("T16189C", 100.0, 22, "mtNG_097")],
                       [mt2("T16189Y", 0.667, "0,1", 16189, "T", "C", "PHP")], weighted=weighted))
    assert list(report.index) == [reported]
    assert report.loc[reported, "Type"] == kind
    assert report.loc[reported, disagreeing] == "DISAGREEMENT"


def test_minor_deletion_from_both_callers_shares_a_row(tmp_path):
    report = calls(run(tmp_path, [fds("A523a", 40.0, 500, "mtNG_005")],
                       [mt2("A523a", 0.40, "300,200", 513, "GCA", "G", "LHP")]))
    assert list(report.index) == ["A523a"]
    assert report.loc["A523a", "called_by_MUTECT2"] is True


def test_deletion_major_minor_split_is_a_disagreement(tmp_path):
    # FDSTOOLS 95% (A523-), Mutect2 85% (A523a): average 90% -> major
    report = calls(run(tmp_path, [fds("A523-", 95.0, 500, "mtNG_005")],
                       [mt2("A523a", 0.85, "75,425", 513, "GCA", "G", "LHP")]))
    assert list(report.index) == ["A523-"]
    assert report.loc["A523-", "Type"] == "DEL"
    assert report.loc["A523-", "called_by_MUTECT2"] == "DISAGREEMENT"


def test_low_rows_per_amplicon_for_both_callers(tmp_path):
    report = run(tmp_path, [fds("T152C", 100.0, 498, "mtNG_001"), fds_low("mtNG_002", 6)], [],
                 depths={"mtNG_003": 3})
    low = report[report["FMP"] == "LOW"].set_index("marker")
    assert sorted(low.index) == ["mtNG_002", "mtNG_003"]
    assert (low.loc["mtNG_002", "called_by_FDSTOOLS"], low.loc["mtNG_002", "called_by_MUTECT2"]) == ("low", "ok")
    assert (low.loc["mtNG_003", "called_by_FDSTOOLS"], low.loc["mtNG_003", "called_by_MUTECT2"]) == ("ok", "low")
    assert low.loc["mtNG_003", "rd_MT2"] == 3


def test_script_writes_excel_with_flags_coloured(tmp_path):
    fds_file, mt2_file, depth_file, library = write_inputs(
        tmp_path,
        [fds("A263G", 100.0, 900, "mtNG_003"), fds("-309.1C", 92.0, 900, "mtNG_003"), fds_low("mtNG_002", 6)],
        [mt2("-309.1c", 0.80, "180,720", 302, "A", "AC", "LHP")], {})
    xlsx = tmp_path / "merged.xlsx"
    subprocess.run([sys.executable, merge.__file__, str(fds_file), str(mt2_file), str(xlsx),
                    "--mutect2_depth", str(depth_file), "--marker_map", str(library)], check=True)
    ws = load_workbook(xlsx).active
    header = [c.value for c in ws[1]]
    col = {name: header.index(name) for name in ("FMP", "called_by_FDSTOOLS", "called_by_MUTECT2")}
    rows = {r[col["FMP"]].value: r for r in ws.iter_rows(min_row=2)}
    fill = lambda cell: cell.fill.start_color.rgb[-6:]
    assert [r[col["FMP"]].value for r in ws.iter_rows(min_row=2)][-1] == "LOW"
    assert fill(rows["A263G"][col["called_by_MUTECT2"]]) == "FFC7CE"          # missed call: red
    assert fill(rows["-309.1c"][col["called_by_FDSTOOLS"]]) == "FFD966"       # disagreement: its own colour
    assert rows["LOW"][col["called_by_FDSTOOLS"]].value == "low"
    assert fill(rows["LOW"][col["called_by_FDSTOOLS"]]) == "FEFE01"          # low coverage: LOW yellow
    assert rows["LOW"][col["called_by_MUTECT2"]].value == "ok"
    assert fill(rows["LOW"][col["called_by_MUTECT2"]]) != "FFC7CE"           # ok is not a miss
