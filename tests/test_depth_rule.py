"""depth_rule.py and resources/noise_model/build_noise_model.py: the depth rule's thresholds from a noise model."""
import importlib.util
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import depth_rule as dr

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "resources" / "fdstools" / "mtNG_library_file.txt"
MODEL = ROOT / "resources" / "noise_model" / "2800M_8controls.tsv"
_spec = importlib.util.spec_from_file_location("build_noise_model", ROOT / "resources" / "noise_model" / "build_noise_model.py")
build = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build)


def model_file(path, rows, masked="8078"):
    """A noise model with one row per control: (name, counts at the 39 shares)."""
    lines = ["# a test model", "\t".join(["control", "median_reads_per_amplicon", "masked"] + dr.SHARE_COLUMNS)]
    lines += ["\t".join([name, "3000", masked] + [str(c) for c in counts]) for name, counts in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def falling(top):
    """Counts that fall from top at 0.1% to 0 at 1.5%."""
    return [max(0, round(top * (0.001 / s) ** 2)) if s < 0.015 else 0 for s in dr.SHARES]


@pytest.mark.parametrize("n", [10, 57, 300])
def test_binomial_tail_is_the_exact_sum(n):
    p = np.array([1e-4, 0.003, 0.02, 0.05])
    for k in range(0, n + 2, max(1, n // 20)):
        exact = [sum(math.comb(n, j) * q ** j * (1 - q) ** (n - j) for j in range(k, n + 1)) for q in p]
        assert dr.binomial_tail(k, n, p) == pytest.approx(exact, rel=1e-9, abs=1e-25)


def test_models_are_pooled(tmp_path, reference):
    a = model_file(tmp_path / "a.tsv", [("c1", falling(300)), ("c2", falling(100))])
    b = model_file(tmp_path / "b.tsv", [("c3", falling(500))], masked="8078,12630")
    counts, masked, names = dr.read_models([a, b])
    assert names == ["c1", "c2", "c3"] and masked == {8078, 12630} and counts.shape == (3, len(dr.SHARES))
    both = model_file(tmp_path / "both.tsv", [("c1", falling(300)), ("c2", falling(100)), ("c3", falling(500))],
                      masked="8078,12630")
    pooled, single = dr.DepthRule([a, b], LIBRARY, reference, 0.01), dr.DepthRule([both], LIBRARY, reference, 0.01)
    assert [pooled.threshold(n) for n in (100, 1000)] == [single.threshold(n) for n in (100, 1000)]


def test_a_file_that_is_not_a_model_is_refused(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("control\tmasked\t0.10\nc1\t\t5\n")
    with pytest.raises(ValueError, match="not a noise model"):
        dr.read_models([bad])


def test_coverage_groups_three_wrong_bases_per_counted_position():
    segments = {"a": [(50, 100)], "b": [(90, 120)]}
    reference = "A" * 16569
    reference = reference[:94] + "N" + reference[95:]                  # 95: no rCRS base
    groups = dr.coverage_groups(segments, reference, masked={110})
    # 50-120 minus the C-stretch 57-60, 95 and the masked 110; 90-100 read by both amplicons
    assert groups == {2: 3 * 10, 1: 3 * (71 - 4 - 1 - 1 - 10)}


def test_shipped_model_thresholds(reference):
    one_in_100 = dr.DepthRule([MODEL], LIBRARY, reference, 0.01)
    one_in_1000 = dr.DepthRule([MODEL], LIBRARY, reference, 0.001)
    depths = [20, 50, 100, 200, 500, 1000, 2000, 5000]
    low, high = [one_in_100.threshold(n) for n in depths], [one_in_1000.threshold(n) for n in depths]
    assert low == SHIPPED_1IN100 and high == SHIPPED_1IN1000
    assert low == sorted(low, reverse=True) and all(h >= l for h, l in zip(high, low))


def sast_file(path, rows):
    pd.DataFrame([{"sample": "ctl", "marker": m, "sequence": s, "flags": "", "total": n} for m, s, n in rows]
                 ).to_csv(path, sep="\t", index=False)
    return path


def test_build_counts_the_wrong_bases_of_a_control(tmp_path, reference):
    segments = {"mtNG_001": [(19, 155)], "mtNG_002": [(133, 267)]}
    sast = sast_file(tmp_path / "ctl.sast.csv", [
        ("mtNG_001", "A73G", 2994),          # the control's own: 73G
        ("mtNG_001", "REF", 6),              # 73A: a wrong base with 6 reads in 3,000
        ("mtNG_001", "A73G C150T", 4),       # 150T: 4 reads in 3,000 + 2,000 (150 is read by both amplicons)
        ("mtNG_001", "A73G T58C", 50),       # in the C-stretch 57-60: not counted
        ("mtNG_002", "REF", 1996),
        ("mtNG_002", "A200G", 4),            # 200G masked: not counted
    ])
    median, counts = build.count_control(sast, segments, "".join(reference), masked={200})
    assert median == 2527                    # of 3,054 and 2,000 reads
    shares = dict(zip(dr.SHARE_COLUMNS, counts))
    # 73A: 6 / 3,054 = 0.196%, 150T: 4 / 5,054 = 0.079%
    assert shares["0.10"] == 1 and shares["0.15"] == 1 and shares["0.20"] == 0


def test_build_writes_a_model_the_rule_reads(tmp_path, reference_fasta):
    sast = sast_file(tmp_path / "ctl.sast.csv", [("mtNG_001", "REF", 2990), ("mtNG_001", "C150T", 10)])
    out = tmp_path / "model.tsv"
    with open(out, "w") as h:
        subprocess.run([sys.executable, str(ROOT / "resources" / "noise_model" / "build_noise_model.py"),
                        "--library", str(LIBRARY), "--reference", str(reference_fasta), str(sast)], stdout=h, check=True)
    counts, masked, names = dr.read_models([out])
    assert names == ["ctl"] and masked == set(build.MASK_2800M) and counts[0][0] == 1


# the shipped model at 20, 50, 100, 200, 500, 1000, 2000 and 5000 reads per amplicon
SHIPPED_1IN100 = [10.0, 8.0, 5.0, 3.5, 2.6, 2.2, 1.95, 1.76]
SHIPPED_1IN1000 = [15.0, 8.0, 6.01, 4.5, 3.2, 2.6, 2.2, 1.92]
