#!/usr/bin/env python3
"""Compare the pipeline's merged calls for simulated samples with their truth.

For every <sample>_truth.tsv (written by simulate_reads.py) the matching
<sample>_merged_variants.xlsx is read. Truth variants are written in the
pipeline's notation (A73G, -315.1C, C8281-; 50/50 heteroplasmies as
T16519Y, -315.1c, C16193c) and compared with the reported
calls (FMP) and with each caller's own label. Calls inside the C-stretches
303-315 and 16180-16193 are counted separately. Mismatches within 10 bp of
each other form one cluster; each cluster is checked on sequence level
(truth vs reported calls applied to rCRS), so "same sequence, different
notation" is told apart from a real error; a heteroplasmy reported as
homoplasmic (or the other way round) counts as a different sequence. --exclude drops truth variants
and calls inside the given windows (e.g. the C-stretches) before comparing.

Writes summary.tsv (one row per sample) and details.tsv (one row per
mismatch) to --out.
"""

import argparse
import re
from pathlib import Path

from openpyxl import load_workbook

REPO = Path(__file__).resolve().parents[2]
C_STRETCHES = [(303, 315), (16180, 16193)]


def load_rcrs(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">"))[:16569].upper()


def expected_labels(variants_used, rcrs):
    labels = set()
    for v in variants_used.split():
        if m := re.fullmatch(r"(\d+)\.(\d+)([ACGTacgt])", v):
            labels.add(f"-{m[1]}.{m[2]}{m[3]}")
        elif m := re.fullmatch(r"(\d+)(-|[ACGTacgtRYMKSW])", v):
            labels.add(f"{rcrs[int(m[1]) - 1]}{m[1]}{m[2]}")
    return labels


def apply_labels(labels, rcrs):
    """Per rCRS position: the bases of the sequence the labels describe.
    Minor forms stay lowercase (A523a, -315.1c) and IUPAC codes stay as
    they are, so a heteroplasmy only matches a heteroplasmy."""
    bases = {p: (b if b in "ACGT" else "") for p, b in enumerate(rcrs, start=1)}
    insertions = {}
    for label in labels:
        label = str(label)
        if m := re.fullmatch(r"-(\d+)\.(\d+)([A-Za-z])", label):
            insertions.setdefault(int(m[1]), {})[int(m[2])] = m[3]
        elif (m := re.fullmatch(r"([ACGTN])(\d+)(-|[A-Za-z])", label)):
            bases[int(m[2])] = "" if m[3] == "-" else m[3]
    for pos, ins in insertions.items():
        bases[pos] += "".join(ins[i] for i in sorted(ins))
    return bases


def clusters(labels, gap=10):
    positions = sorted({position(l) for l in labels if position(l)})
    groups = []
    for p in positions:
        if groups and p - groups[-1][-1] <= gap:
            groups[-1].append(p)
        else:
            groups.append([p])
    return [(g[0], g[-1]) for g in groups]


def position(label):
    m = re.search(r"(\d+)", str(label))
    return int(m[1]) if m else None


def region(label):
    p = position(label)
    return next((f"{a}-{b}" for a, b in C_STRETCHES if p and a <= p <= b), "elsewhere")


def read_report(xlsx):
    ws = load_workbook(xlsx).active
    header = [c.value for c in ws[1]]
    if "FMP" not in header:
        return [], []
    col = {h: i for i, h in enumerate(header)}
    calls, low = [], []
    for row in ws.iter_rows(min_row=2, values_only=True):
        get = lambda k: row[col[k]] if k in col else None
        if get("FMP") == "LOW":
            low.append(f"{get('marker')}(FDS={get('called_by_FDSTOOLS')},MT2={get('called_by_MUTECT2')})")
        else:
            calls.append({k: get(k) for k in ("FMP", "FDSTOOLS", "MUTECT2", "called_by_FDSTOOLS", "called_by_MUTECT2")})
    return calls, low


def compare(expected, labels, rcrs):
    """Missed and extra labels, and per mismatch cluster whether the sequence
    these labels describe equals the true sequence."""
    missed, extra = expected - labels, labels - expected
    true_seq, seq = apply_labels(expected, rcrs), apply_labels(labels, rcrs)
    same = {}
    for start, end in clusters(missed | extra):
        window = range(max(1, start - 5), min(len(rcrs), end + 5) + 1)
        ok = "".join(true_seq[p] for p in window) == "".join(seq[p] for p in window)
        for p in range(start, end + 1):
            same[p] = ok
    wrong = sum(1 for c in clusters(missed | extra) if not same[c[0]])
    return missed, extra, same, wrong


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--merged-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--reference", default=REPO / "resources/rCRS/rCRS_NimaGen.fasta")
    p.add_argument("--exclude", default="", help="windows to leave out, e.g. 300-320,16170-16200")
    args = p.parse_args()
    windows = [tuple(map(int, w.split("-"))) for w in args.exclude.split(",") if w]
    keep = lambda labels: {l for l in labels if not any(a <= (position(l) or 0) <= b for a, b in windows)}

    rcrs = load_rcrs(args.reference)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    views = ("report", "fdstools", "mutect2")
    summary, details = [], []
    for truth_file in sorted(Path(args.truth_dir).glob("*_truth.tsv")):
        sample = truth_file.name[: -len("_truth.tsv")]
        truth = dict(l.rstrip("\n").split("\t", 1) for l in open(truth_file))
        expected = keep(expected_labels(truth.get("variants_used", ""), rcrs))
        xlsx = Path(args.merged_dir) / f"{sample}_merged_variants.xlsx"
        if not xlsx.exists():
            summary.append([sample, truth["haplogroup"], len(expected), "NO_REPORT"])
            continue
        calls, low = read_report(xlsx)
        labels = {
            "report": keep({c["FMP"] for c in calls}),
            "fdstools": keep({c["FDSTOOLS"] for c in calls if c["FDSTOOLS"]}),
            "mutect2": keep({c["MUTECT2"] for c in calls if c["MUTECT2"]}),
        }
        row = [sample, truth["haplogroup"], len(expected)]
        for view in views:
            missed, extra, same, wrong = compare(expected, labels[view], rcrs)
            row += [len(missed), len(extra), wrong]
            for kind, group in (("missed", missed), ("extra", extra)):
                for label in sorted(group, key=position):
                    details.append([sample, truth["haplogroup"], view, kind, label, region(label),
                                    "same sequence" if same.get(position(label)) else "sequence differs"])
        summary.append(row + [" ".join(low)])

    with open(out / "summary.tsv", "w") as f:
        cols = [f"{v}_{k}" for v in views for k in ("missed", "extra", "wrong_sequence_clusters")]
        f.write("\t".join(["sample", "haplogroup", "expected"] + cols + ["low_amplicons"]) + "\n")
        for row in summary:
            f.write("\t".join(map(str, row)) + "\n")
    with open(out / "details.tsv", "w") as f:
        f.write("sample\thaplogroup\tview\tkind\tlabel\tregion\tsequence_check\n")
        for row in details:
            f.write("\t".join(map(str, row)) + "\n")
    scored = [r for r in summary if r[3] != "NO_REPORT"]
    for i, view in enumerate(views):
        exact = sum(1 for r in scored if r[3 + 3 * i] == 0 and r[4 + 3 * i] == 0)
        right_seq = sum(1 for r in scored if r[5 + 3 * i] == 0)
        print(f"{view:9}: {exact}/{len(scored)} samples exact, {right_seq}/{len(scored)} with the correct sequence")


if __name__ == "__main__":
    main()
