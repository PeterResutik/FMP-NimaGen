"""merge_fdstools_mutect2_improved.py: FDSTOOLS' table -> the report."""
import subprocess
import sys

import pandas as pd
import pytest
from openpyxl import load_workbook

import merge_fdstools_mutect2_improved as merge

AMPLICONS = {"mtNG_001": (19, 155), "mtNG_002": (133, 267), "mtNG_003": (259, 367),
             "mtNG_005": (470, 610), "mtNG_097": (16094, 16276)}


def fds(label, vf, depth, marker):
    """A row of the FDSTOOLS table (vf in percent)."""
    start, end = AMPLICONS[marker]
    return {"FDSTOOLS": label, "vf_FDS": vf, "rd_FDS": round(vf * depth / 100),
            "rd_total_FDS": depth, "marker": marker,
            "marker_range": f"{start}-{end}", "num_markers": 1}


def fds_bases(label, vf, reads, depth, marker):
    """An FDSTOOLS row whose vf and reads are given per base ("C 40, A 10")."""
    return {**fds(label, 0.0, depth, marker), "vf_FDS": vf, "rd_FDS": reads}


def fds_low(marker, reads):
    start, end = AMPLICONS[marker]
    return {"FDSTOOLS": "LOW", "rd_FDS": reads, "marker": marker, "marker_range": f"{start}-{end}"}


def write_inputs(tmp_path, fds_rows):
    library = tmp_path / "library.txt"
    library.write_text("[genome_position]\n" + "".join(
        f"{m} = chrM, {s}, {e}\n" for m, (s, e) in AMPLICONS.items()) + "[flanks]\n")
    fds_file = tmp_path / "fds.txt"
    pd.DataFrame(fds_rows).to_csv(fds_file, sep="\t", index=False)
    return fds_file, library


def run(tmp_path, fds_rows, min_vf=5.0):
    fds_file, library = write_inputs(tmp_path, fds_rows)
    return merge.fdstools_report(str(fds_file), min_vf=min_vf, marker_map=str(library))


def calls(report):
    return report[report["FMP"] != "LOW"].set_index("FMP")


def run_script(tmp_path, fds_rows, *options):
    fds_file, library = write_inputs(tmp_path, fds_rows)
    xlsx = tmp_path / "report.xlsx"
    subprocess.run([sys.executable, merge.__file__, str(fds_file), str(xlsx), "--marker_map", str(library),
                    *options], check=True)
    ws = load_workbook(xlsx).active
    return [c.value for c in ws[1]], list(ws.iter_rows(min_row=2))


@pytest.mark.parametrize("a, b", [("T16189C", "T16189Y"), ("T16189Y", "T16189H"), ("C756M", "C756Y"),
                                  ("A523-", "A523a")])
def test_order_key_same_locus_any_spelling(a, b):
    assert merge.order_key(a) == merge.order_key(b)


def test_rows_by_position_deletion_before_substitution(tmp_path):
    report = calls(run(tmp_path, [fds("A16183C", 100.0, 500, "mtNG_097"), fds("A263G", 100.0, 900, "mtNG_003"),
                                  fds("-309.1C", 92.0, 900, "mtNG_003"), fds("A16183-", 40.0, 500, "mtNG_097")]))
    assert list(report.index) == ["A263G", "-309.1C", "A16183-", "A16183C"]


def test_report_holds_fdstools_columns_only(tmp_path):
    # the columns of the FDSTOOLS step's table, in its order
    columns = ["FDSTOOLS", "vf_FDS", "rd_FDS", "rd_total_FDS", "variant_note", "dominant_molecule",
               "marker", "marker_range", "num_markers", "is_noise_or_low_frq"]
    row = {**fds("T152C", 100.0, 498, "mtNG_001"), "variant_note": None, "dominant_molecule": None,
           "is_noise_or_low_frq": False}
    report = run(tmp_path, [{c: row[c] for c in columns}])
    assert list(report.columns) == ["FMP", "vf_FDS", "rd_FDS", "rd_total_FDS", "variant_note", "dominant_molecule",
                                    "marker", "marker_range", "num_markers", "is_noise_or_low_frq"]


def test_rows_per_base_become_one_row(tmp_path):
    # T 50%, C 40%, A 10% written as one row per base other than rCRS
    report = calls(run(tmp_path, [fds("T16189Y", 40.0, 500, "mtNG_097"), fds("T16189W", 10.0, 500, "mtNG_097")]))
    assert list(report.index) == ["T16189H"]
    assert (report.loc["T16189H", "vf_FDS"], report.loc["T16189H", "rd_FDS"]) == ("C 40, A 10", "C 200, A 50")


def test_rows_per_base_without_rcrs(tmp_path):
    # C 70% and A 30%, no T: the old rows' codes both claim a T
    report = calls(run(tmp_path, [fds("T16189Y", 70.0, 500, "mtNG_097"), fds("T16189W", 30.0, 500, "mtNG_097")]))
    assert list(report.index) == ["T16189M"]
    assert report.loc["T16189M", "vf_FDS"] == "C 70, A 30"


@pytest.mark.parametrize("deleted, reported", [(None, "T16189H"), (10.0, "T16189M")])
def test_rcrs_share_leaves_out_deleted_molecules(tmp_path, deleted, reported):
    # C 50%, A 40%: the remaining 10% is rCRS T, or molecules without a base at 16189
    rows = [fds("T16189Y", 50.0, 500, "mtNG_097"), fds("T16189W", 40.0, 500, "mtNG_097")]
    if deleted:
        rows.append(fds("T16189t", deleted, 500, "mtNG_097"))
    assert reported in calls(run(tmp_path, rows)).index


def test_rows_given_per_base_are_left_as_they_are(tmp_path):
    report = calls(run(tmp_path, [fds_bases("A16265V", "C 6, G 4", "C 30, G 20", 500, "mtNG_097"),
                                  fds("A16265G", 3.0, 500, "mtNG_097"), fds("-309.1C", 92.0, 900, "mtNG_003")]))
    assert {"A16265V", "A16265G", "-309.1C"} <= set(report.index)
    assert report.loc["A16265V", "vf_FDS"] == "C 6, G 4"


def test_dominant_molecule_passes_through(tmp_path):
    report = calls(run(tmp_path, [{**fds("T16189C", 100.0, 900, "mtNG_097"), "dominant_molecule": "T16189C (62.0%)"}]))
    assert report.loc["T16189C", "dominant_molecule"] == "T16189C (62.0%)"


def test_low_rows_per_amplicon_at_the_end(tmp_path):
    report = run(tmp_path, [fds_low("mtNG_097", 4), fds("T152C", 100.0, 498, "mtNG_001"), fds_low("mtNG_002", 6)])
    low = report[report["FMP"] == "LOW"]
    assert list(report["FMP"]) == ["T152C", "LOW", "LOW"]
    assert list(low["marker"]) == ["mtNG_002", "mtNG_097"]      # in the library's order
    assert list(low["rd_FDS"]) == [6, 4]
    assert list(low["marker_range"]) == ["133-267", "16094-16276"]


def test_empty_table_gives_an_empty_report(tmp_path):
    fds_file, library = write_inputs(tmp_path, [])
    fds_file.write_text("FDSTOOLS\tvf_FDS\trd_FDS\trd_total_FDS\tmarker\tmarker_range\tnum_markers\n")
    assert merge.fdstools_report(str(fds_file), marker_map=str(library)).empty


def test_script_writes_excel_coloured_by_kind(tmp_path):
    header, rows = run_script(tmp_path, [fds("A263G", 100.0, 900, "mtNG_003"), fds("-309.1c", 30.0, 900, "mtNG_003"),
                                         fds_bases("T16189H", "C 40, A 10", "C 200, A 50", 500, "mtNG_097"),
                                         fds_low("mtNG_002", 6)])
    fill = lambda cell: cell.fill.start_color.rgb[-6:]
    by_call = {r[0].value: r for r in rows}
    assert "MUTECT2" not in header and "called_by_FDSTOOLS" not in header and "FDSTOOLS" not in header
    assert rows[-1][0].value == "LOW" and fill(rows[-1][0]) == "FEFE01"      # LOW: yellow, last
    assert fill(by_call["T16189H"][0]) == "FAC000"                             # an IUPAC code
    assert fill(by_call["-309.1c"][0]) == "3CB0F1"                             # lowercase: a length variant
    assert fill(by_call["A263G"][header.index("vf_FDS")]) == "C6EFCE"          # FDSTOOLS' columns: green
    assert by_call["T16189H"][header.index("vf_FDS")].value == "C 40, A 10"


@pytest.mark.parametrize("variant_note, written", [("true", True), ("false", False)])
def test_script_leaves_out_the_variant_note(tmp_path, variant_note, written):
    header, _rows = run_script(tmp_path, [{**fds("T16189C", 100.0, 900, "mtNG_097"), "variant_note": "a note"}],
                               "--variant_note", variant_note)
    assert ("variant_note" in header) == written
