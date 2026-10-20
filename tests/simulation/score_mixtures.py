#!/usr/bin/env python3
"""Score the reports of simulated two-haplogroup mixtures (simulate_reads.py --mixtures).

Sites are the positions where the two haplogroups carry different single bases, outside the C-stretches 303-315
and 16180-16193 and more than 5 bases away from an insertion or deletion of either haplogroup. At every site the
minor haplogroup's base is on the given share of each amplicon's reads. Per site and per view (the report's FMP
label, FDSTOOLS' and Mutect2's own labels):
  minor found  the label is the IUPAC code of both bases (T1119Y);
  major only   the label is the major base, or there is none where the major base is rCRS's;
  other        any other label.
The expected frequency of a site is that of its non-rCRS base(s), as the report gives it (vf_FDS, vf_MT2).
Calls at positions where neither haplogroup differs from rCRS (C-stretches aside) are extra calls.

Writes summary.tsv (one row per mixture) and sites.tsv (one row per site) to --out.
"""
import argparse
import re
from pathlib import Path
from statistics import median

from openpyxl import load_workbook

from score_calls import REPO, C_STRETCHES, apply_labels, expected_labels, load_rcrs, position

IUPAC_CODE = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("AC"): "M", frozenset("GT"): "K",
              frozenset("CG"): "S", frozenset("AT"): "W"}
SUBSTITUTION = re.compile(r"[ACGTN]\d+[ACGTRYMKSW]")
VIEWS = (("report", "FMP", None), ("fdstools", "FDSTOOLS", "vf_FDS"), ("mutect2", "MUTECT2", "vf_MT2"))


def in_c_stretch(p):
    return any(a <= p <= b for a, b in C_STRETCHES)


def sequence(variants_used, rcrs):
    """Per rCRS position the haplogroup's base(s), from the truth file's variants."""
    return apply_labels(expected_labels(variants_used, rcrs), rcrs)


def sites(major, minor, rcrs):
    """[(position, major base, minor base)] where the two differ by one base each, away from indels."""
    indels = {p for seq in (major, minor) for p, b in seq.items() if len(b) != 1}
    near_indel = {q for p in indels for q in range(p - 5, p + 6)}
    return [(p, major[p], minor[p]) for p in sorted(major)
            if major[p] != minor[p] and len(major[p]) == 1 and len(minor[p]) == 1
            and major[p] in "ACGT" and minor[p] in "ACGT" and not in_c_stretch(p) and p not in near_indel]


def status(label, ref, p, base_major, base_minor):
    if label == f"{ref}{p}{IUPAC_CODE[frozenset((base_major, base_minor))]}":
        return "minor found"
    if label == (f"{ref}{p}{base_major}" if base_major != ref else None):
        return "major only"
    return "other"


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_rows(xlsx):
    ws = load_workbook(xlsx).active
    header = [c.value for c in ws[1]]
    return [dict(zip(header, row)) for row in ws.iter_rows(min_row=2, values_only=True)
            if row and row[0] is not None and row[header.index("FMP")] != "LOW"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--merged-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--reference", default=REPO / "resources/rCRS/rCRS_NimaGen.fasta")
    p.add_argument("--frame", choices=["separate", "shared"], default="separate")
    args = p.parse_args()

    rcrs = load_rcrs(args.reference)
    plain = apply_labels(set(), rcrs)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    summary, site_rows = [], []
    for truth_file in sorted(Path(args.truth_dir).glob("*_truth.tsv")):
        truth = dict(l.rstrip("\n").split("\t", 1) for l in open(truth_file))
        if "minor_haplogroup" not in truth:
            continue
        sample = truth_file.name[: -len("_truth.tsv")]
        head = [sample, truth["haplogroup"], truth["minor_haplogroup"], truth["minor_share"],
                truth.get("reads_per_amplicon", ""), truth.get("error_source", "")]
        xlsx = Path(args.merged_dir) / f"{sample}_{args.frame}_frame.xlsx"
        if not xlsx.exists():
            summary.append(head + ["NO_REPORT"])
            continue
        major = sequence(truth.get("variants_used", ""), rcrs)
        minor = sequence(truth.get("minor_variants_used", ""), rcrs)
        share = float(truth["minor_share"]) / 100
        rows = read_rows(xlsx)
        by_pos = {}
        for r in rows:
            if SUBSTITUTION.fullmatch(str(r["FMP"])):
                by_pos.setdefault(position(r["FMP"]), r)
        found = {v: 0 for v, _, _ in VIEWS}
        counts = {"major only": 0, "other": 0}
        errors = {"fdstools": [], "mutect2": []}
        site_list = sites(major, minor, rcrs)
        for pos, base_major, base_minor in site_list:
            ref = rcrs[pos - 1]
            r = by_pos.get(pos, {})
            expected_vf = 100 * ((base_major != ref) * (1 - share) + (base_minor != ref) * share)
            row = [sample, truth["haplogroup"], truth["minor_haplogroup"], truth["minor_share"],
                   truth.get("reads_per_amplicon", ""), truth.get("error_source", ""), pos, ref, base_major,
                   base_minor, round(expected_vf, 2)]
            for view, col, vf_col in VIEWS:
                st = status(r.get(col), ref, pos, base_major, base_minor)
                row += [r.get(col) or "", st]
                if st == "minor found":
                    found[view] += 1
                    vf = number(r.get(vf_col)) if vf_col else None
                    if vf is not None:
                        errors[view].append(vf - expected_vf)
                elif view == "report":
                    counts[st] += 1
                if vf_col:
                    row.append(r.get(vf_col) if r.get(vf_col) is not None else "")
            site_rows.append(row)
        explained = {q for q, b in major.items() if b != plain[q]} | {q for q, b in minor.items() if b != plain[q]}
        extra = sum(1 for r in rows if position(r["FMP"]) is not None and not in_c_stretch(position(r["FMP"]))
                    and position(r["FMP"]) not in explained)
        stats = []
        for view in ("fdstools", "mutect2"):
            e = errors[view]
            stats += [round(median(abs(x) for x in e), 2) if e else "", round(median(e), 2) if e else ""]
        summary.append(head + [len(site_list), found["report"], found["fdstools"], found["mutect2"],
                               counts["major only"], counts["other"], extra] + stats)

    with open(out / "summary.tsv", "w") as f:
        f.write("\t".join(["sample", "major", "minor", "minor_share_pct", "reads_per_amplicon", "error_source",
                           "sites", "minor_found", "minor_found_FDS", "minor_found_MT2", "major_only", "other",
                           "extra_calls", "abs_err_vf_FDS", "bias_vf_FDS", "abs_err_vf_MT2", "bias_vf_MT2"]) + "\n")
        for row in summary:
            f.write("\t".join(map(str, row)) + "\n")
    with open(out / "sites.tsv", "w") as f:
        f.write("\t".join(["sample", "major", "minor", "minor_share_pct", "reads_per_amplicon", "error_source",
                           "pos", "rCRS", "base_major", "base_minor", "expected_vf_pct", "FMP", "report_status",
                           "FDSTOOLS", "fdstools_status", "vf_FDS", "MUTECT2", "mutect2_status", "vf_MT2"]) + "\n")
        for row in site_rows:
            f.write("\t".join(map(str, row)) + "\n")
    scored = [r for r in summary if r[6] != "NO_REPORT"]
    total = sum(r[6] for r in scored)
    print(f"{len(scored)} mixtures, {total} sites; minor base found by the report at "
          f"{sum(r[7] for r in scored)}, FDSTOOLS {sum(r[8] for r in scored)}, Mutect2 {sum(r[9] for r in scored)}")


if __name__ == "__main__":
    main()
