"""process_mutect2_output_improved.py: Mutect2 records -> report labels."""
import subprocess
from pathlib import Path
import sys

import pandas as pd
import pytest

import process_mutect2_output_improved as mt2


@pytest.mark.parametrize("a, b", [
    ("T16189C", "T16189Y"),
    ("-309.1C", "-309.1c"),
    ("A523-", "A523a"),
])
def test_claim_ignores_major_minor_spelling(a, b):
    assert mt2.claim_key(a) == mt2.claim_key(b)


def test_claim_keeps_different_alleles_apart():
    assert mt2.claim_key("C756M") != mt2.claim_key("C756Y")


@pytest.mark.parametrize("pos, expected", [(16569, 16569), (16570, 1), (16579, 10)])
def test_origin_copy_positions_wrap(pos, expected):
    assert mt2.wrap_circular_position(pos) == expected


def test_snp_in_origin_copy_labelled_at_true_position(reference):
    assert mt2.apply_snp(16579, "T", "C", 1.0, reference, 0.05) == ("T10C", "SNP")


def test_origin_overlap_pooled_like_one_pileup():
    # chrM:21 seen through the appended copy (16590) and at its true position
    df = pd.DataFrame({
        "Pos": [16590, 21], "Ref": ["A", "A"], "Variant": ["G", "G"],
        "Coverage": ["10,90", "20,80"], "VariantLevel": [0.9, 0.8],
    })
    pooled = mt2.pool_origin_overlap(df)
    assert len(pooled) == 1
    assert pooled.loc[0, "Pos"] == 21
    assert pooled.loc[0, "Coverage"] == "30,170"
    assert pooled.loc[0, "VariantLevel"] == pytest.approx(0.85)


@pytest.mark.parametrize("threshold", [0.10, 0.90])
def test_lh_bounds_same_for_either_side(threshold):
    assert mt2.lh_bounds(threshold) == pytest.approx((0.10, 0.90))


@pytest.mark.parametrize("level, expected", [
    (0.05, (None, "BELOW_LH_FLOOR")),
    (0.50, ("-309.1c", "LHP")),
    (0.95, ("-309.1C", "INS")),
])
def test_insertion_floor_minor_major(reference, level, expected):
    # a C inserted in the 303-309 C-run, left-aligned by Mutect2 after A302
    assert mt2.apply_insertion(302, "A", "AC", level, reference, 0.10) == expected


@pytest.mark.parametrize("level, expected", [
    (0.05, (None, "BELOW_LH_FLOOR")),
    (0.50, ("A523a C524c", "LHP")),
    (1.00, ("A523- C524-", "DEL")),
])
def test_deletion_floor_minor_major(reference, level, expected):
    # one AC unit of the 514-524 repeat, left-aligned by Mutect2 after G513
    assert mt2.apply_deletion(513, "GCA", "G", level, reference, 0.10) == expected


def test_repeat_deletion_lands_on_last_copy(reference):
    # CCCCCTCTA deleted, left-aligned at 8271: belongs at 8281-8289, not 8280-8288
    assert mt2.apply_deletion(8271, "ACCCCCTCTA", "A", 1.0, reference, 0.10) == (
        "C8281- C8282- C8283- C8284- C8285- T8286- C8287- T8288- A8289-", "DEL")


@pytest.mark.parametrize("level, expected", [(0.293, ("C756M", "PHP")), (0.97, ("C756A", "SNP"))])
def test_substitution_minor_iupac_major_base(reference, level, expected):
    assert mt2.apply_snp(756, "C", "A", level, reference, 0.05) == expected


def test_script_end_to_end(tmp_path, reference_fasta):
    """The table p13 builds from the VCF, through the whole script."""
    columns = ["ID", "Filter", "Pos", "Ref", "Variant", "VariantLevel", "MeanBaseQuality", "Coverage", "GT", "Type"]
    rows = [
        ["s.bam", "PASS", 21, "A", "G", 0.9, "35,35", "10,90", "0/1", "SNP"],
        ["s.bam", "PASS", 756, "C", "A", 0.293, "35,35", "70,29", "0/1", "SNP"],
        ["s.bam", "PASS", 756, "C", "T", 0.051, "35,35", "94,5", "0/1", "SNP"],
        ["s.bam", "PASS", 16590, "A", "G", 0.8, "35,35", "20,80", "0/1", "SNP"],
    ]
    table = tmp_path / "mutect2.txt"
    pd.DataFrame(rows, columns=columns).to_csv(table, sep="\t", index=False)
    out = tmp_path / "mutect2.empop.txt"
    subprocess.run([sys.executable, mt2.__file__, str(table), str(out), str(reference_fasta),
                    "--min_vf", "5", "--lh_thresh", "10"], check=True)
    result = pd.read_csv(out, sep="\t")
    # two bases other than rCRS at 756 become one row; chrM:21 seen twice one pooled row
    calls = dict(zip(result["MUTECT2"], result["vf_MT2"]))
    assert sorted(calls) == ["A21R", "C756H"]
    assert float(calls["A21R"]) == pytest.approx(0.85)
    assert calls["C756H"] == "A 0.293, T 0.051"
    assert result.set_index("MUTECT2").loc["A21R", "rd_MT2"] == "30,170"
    assert result.set_index("MUTECT2").loc["C756H", "rd_MT2"] == "70,29,5"
    assert result.set_index("MUTECT2").loc["C756H", "Type"] == "PHP"


def test_script_writes_major_calls_by_the_general_rule(tmp_path, reference_fasta):
    """L3f3: Mutect2 records one T deleted and T15941C; the same molecule is
    T15940C T15944- by the general rule, as FDSTOOLS' spelling becomes too."""
    columns = ["ID", "Filter", "Pos", "Ref", "Variant", "VariantLevel", "MeanBaseQuality", "Coverage", "GT", "Type"]
    rows = [
        ["s.bam", "PASS", 15939, "CT", "C", 1.0, "35,35", "0,300", "1/1", "INDEL"],
        ["s.bam", "PASS", 15941, "T", "C", 1.0, "35,35", "0,300", "1/1", "SNP"],
    ]
    table = tmp_path / "mutect2.txt"
    pd.DataFrame(rows, columns=columns).to_csv(table, sep="\t", index=False)
    out = tmp_path / "mutect2.empop.txt"
    subprocess.run([sys.executable, mt2.__file__, str(table), str(out), str(reference_fasta),
                    "--min_vf", "5", "--lh_thresh", "10"], check=True)
    result = pd.read_csv(out, sep="\t").set_index("MUTECT2")
    assert sorted(result.index) == ["T15940C", "T15944-"]
    assert (result.loc["T15940C", "Type"], result.loc["T15944-", "Type"]) == ("SNP", "DEL")


def test_a_label_the_rewrite_keeps_keeps_its_frequency(reference):
    # G15933A stays as it is; only T15940- T15941C become T15940C T15944-
    df = pd.DataFrame([{"MUTECT2": l, "VariantLevel": v, "Coverage": "0,100", "Type": "?"}
                       for l, v in (("G15933A", 1.0), ("T15940-", 0.96), ("T15941C", 0.97))])
    out = mt2.respell_rows(df, reference).set_index("MUTECT2")["VariantLevel"].to_dict()
    assert out == {"G15933A": 1.0, "T15940C": 0.96, "T15944-": 0.96}


def frame(labels, which):
    df = pd.DataFrame([{"MUTECT2": l, "VariantLevel": v, "Coverage": "0,100", "Type": "?"} for l, v in labels])
    return mt2.frame_rows(df, mt2.load_reference(REFERENCE), which).set_index("MUTECT2")["VariantLevel"].to_dict()


REFERENCE = str(Path(__file__).resolve().parents[1] / "resources" / "rCRS" / "rCRS_NimaGen.fasta")


@pytest.mark.parametrize("labels, which, expected", [
    # D5c as Mutect2 records it; the shared frame writes it as FDSTOOLS' molecules do
    ([("A16181-", .99), ("A16182-", .99), ("A16183-", .99), ("T16189C", .99), ("-16192.1T", .99),
      ("-16192.2C", .99), ("-16192.3C", .99)], "shared",
     {"A16181C": .99, "A16182C": .99, "A16183C": .99, "T16189C": .99, "C16190T": .99}),
    # T57C with one T more (and 55.1T): -56.1C becomes T57C ... -60.1T
    ([("-55.1T", .99), ("-56.1C", .99), ("T58C", .99)], "separate",
     {"-55.1T": .99, "T57C": .99, "T59C": .99, "-60.1T": .99}),
    # next to 57-60 the edges follow the general rule, as for FDSTOOLS
    ([("-60.1T", .99), ("-60.2T", .99), ("C64-", .99), ("T65-", .99)], "shared",
     {"C61T": .99, "G62-": .99, "-64.1G": .99}),
    # a boundary C below the substitution line (95%) is minor and stays as Mutect2 wrote it
    ([("A16183M", .92), ("T16189C", 1.0)], "separate", {"A16183M": .92, "T16189C": 1.0}),
    # a major one is written in the frame; T16189C keeps its own frequency
    ([("A16183C", .97), ("T16189C", 1.0)], "separate", {"A16183-": .97, "T16189C": 1.0, "-16193.1C": .97}),
    ([("A16183M", .86), ("T16189C", 1.0), ("-16193.1c", .54)], "separate",
     {"A16183M": .86, "T16189C": 1.0, "-16193.1c": .54}),
])
def test_mutect2_major_calls_in_the_frames(labels, which, expected):
    assert frame(labels, which) == pytest.approx(expected)


def records(*rows):
    """Mutect2's table after finalize_output_table: (label, frequency, ref,alt reads)."""
    return pd.DataFrame([{"MUTECT2": l, "VariantLevel": v, "Coverage": ad, "Type": "?", "Pos": 0}
                         for l, v, ad in rows])


@pytest.mark.parametrize("rows, expected", [
    # A and G, no C left: the code of A and G only
    ([("C756M", 0.7, "0,70"), ("C756S", 0.3, "0,30")], {"C756R": "A 0.7, G 0.3"}),
    # A on 94%, rCRS on 3 of 97 reads, and a deletion: major
    ([("C756M", 0.94, "3,94"), ("C756c", 0.03, "97,3")], {"C756A": 0.94, "C756c": 0.03}),
    # one base other than rCRS stays as it is
    ([("C756M", 0.293, "70,29")], {"C756M": 0.293}),
])
def test_one_row_per_position(rows, expected):
    result = mt2.position_rows(records(*rows), 0.05)
    assert dict(zip(result["MUTECT2"], result["VariantLevel"])) == expected


@pytest.mark.parametrize("rows, expected", [
    # 100 reads, all G: Mutect2 reports 0.99; its reads show no rCRS A
    ([("A73R", 0.99, "0,100")], {"A73G": 0.99}),
    # 50% rCRS, 30% T and 20% deleted (the synthetic test): Mutect2 counts the deleted
    # reads toward T (0.418), but its reads still give rCRS 300 of 600
    ([("A7025W", 0.418, "300,300"), ("A7025a", 0.201, "480,120")], {"A7025W": 0.418, "A7025a": 0.201}),
])
def test_rcrs_share_from_the_reads(rows, expected):
    result = mt2.position_rows(records(*rows), 0.01)
    assert dict(zip(result["MUTECT2"], result["VariantLevel"])) == expected
