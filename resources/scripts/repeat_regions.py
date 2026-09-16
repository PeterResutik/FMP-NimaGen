#!/usr/bin/env python3
"""Reference-derived repeat regions.

Identifies stretches of the reference where indel placement is
inherently ambiguous - homopolymer runs, and groups of runs close enough
that a single base change can merge them (e.g. rCRS 16184-16193,
CCCCC-T-CCCC, which becomes one 10bp C-run the moment T16189C is
present, as it is in many samples).

These regions are a property of the REFERENCE alone - deterministic,
sample-independent, computed once. That matters: the failure mode this
replaces was hardcoding coordinates inferred from one sample's evidence.
"""

import argparse


def find_homopolymer_runs(seq, min_length=4):
    """Maximal runs of one base, as 1-based (start, end, base, length)."""
    runs = []
    i = 0
    n = len(seq)
    while i < n:
        j = i
        while j + 1 < n and seq[j + 1] == seq[i]:
            j += 1
        length = j - i + 1
        if length >= min_length and seq[i] in "ACGT":
            runs.append((i + 1, j + 1, seq[i], length))
        i = j + 1
    return runs


def merge_runs(runs, max_gap=1):
    """Merge runs separated by <= max_gap reference bases into regions.

    A single interrupting base (max_gap=1) is the important case: a lone
    mismatch between two runs of the same base disappears whenever that
    position is substituted, joining them into one longer run, and indel
    placement then becomes ambiguous across the whole span.
    """
    if not runs:
        return []
    regions = []
    cur_start, cur_end = runs[0][0], runs[0][1]
    members = [runs[0]]
    for start, end, base, length in runs[1:]:
        if start - cur_end - 1 <= max_gap:
            cur_end = end
            members.append((start, end, base, length))
        else:
            regions.append((cur_start, cur_end, members))
            cur_start, cur_end = start, end
            members = [(start, end, base, length)]
    regions.append((cur_start, cur_end, members))
    return regions


def find_boundary_run_regions(seq, min_length=4, max_gap=1):
    """Reference-derived leading-run/extension-run boundary regions - a
    homopolymer run immediately followed by a run of a DIFFERENT base
    (e.g. rCRS 16180-16193, the A-run immediately followed by the
    T-interrupted C-run). Same reference-only, sample-independent
    determinism as find_homopolymer_runs/merge_runs above.

    Restricted to the exact shape validated against real data (see
    report_boundary_run in process_fdstools_output_improved_better.py):
    the LEADING run itself must be one single gapless homopolymer run (no
    internal base-identity change), immediately adjacent to an extension
    run that is entirely ONE other base (a bridging single-base gap
    WITHIN the extension run, e.g. 16184-16193's own T16189, is fine -
    that's the case this was built for). A merged region not matching
    this exactly (a gap between the leading and extension runs, a base
    change inside the leading run itself, or more than two distinct
    bases overall) is left alone entirely, for the general per-sample
    reconciliation to handle as it always has - not force-fit into a
    shape it wasn't checked against.
    """
    regions = merge_runs(find_homopolymer_runs(seq, min_length), max_gap)
    boundary_regions = []
    for start, end, members in regions:
        if len(members) < 2:
            continue
        m0, m1 = members[0], members[1]
        if m0[2] == m1[2]:
            continue
        if m0[1] + 1 != m1[0]:
            continue
        later_bases = {m[2] for m in members[1:]}
        if len(later_bases) != 1:
            continue
        boundary_regions.append({"leading": (m0[0], m0[1]), "extension": (m0[1] + 1, end)})
    return boundary_regions


def load_reference(fasta_path):
    with open(fasta_path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", help="Reference FASTA")
    parser.add_argument("--min-run", type=int, default=4)
    parser.add_argument("--max-gap", type=int, default=1)
    parser.add_argument("--show", action="store_true", help="List every region")
    args = parser.parse_args()

    seq = load_reference(args.reference).upper()
    runs = find_homopolymer_runs(seq, args.min_run)
    regions = merge_runs(runs, args.max_gap)

    multi = [r for r in regions if len(r[2]) > 1]
    print(f"reference length: {len(seq)}")
    print(f"homopolymer runs >= {args.min_run}bp: {len(runs)}")
    print(f"merged regions (gap <= {args.max_gap}): {len(regions)}  "
          f"({len(multi)} of them span more than one run)")
    total_bp = sum(end - start + 1 for start, end, _ in regions)
    print(f"total bases covered: {total_bp} ({100.0*total_bp/len(seq):.2f}% of reference)")

    if args.show:
        for start, end, members in regions:
            desc = " ".join(f"{b}x{l}@{s}" for s, e, b, l in members)
            print(f"  {start}-{end}  ({end-start+1}bp)  {desc}")


if __name__ == "__main__":
    main()
