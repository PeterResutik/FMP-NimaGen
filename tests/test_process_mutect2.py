"""process_mutect2_output_improved.py: Mutect2 records -> report labels."""
import subprocess
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
    # two alleles at 756 stay two rows; chrM:21 seen twice becomes one pooled row
    assert dict(zip(result["MUTECT2"], result["vf_MT2"])) == pytest.approx(
        {"A21R": 0.85, "C756M": 0.293, "C756Y": 0.051})
    assert result.set_index("MUTECT2").loc["A21R", "rd_MT2"] == "30,170"


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
