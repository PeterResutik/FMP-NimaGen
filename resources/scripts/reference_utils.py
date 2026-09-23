#!/usr/bin/env python3
"""Reference-sequence math - everything computed purely from the
reference FASTA itself, deterministic and sample-independent.

Two related groups of primitives:
  - Repeat-region detection: homopolymer runs, merged regions, and
    boundary-run regions (a leading run immediately followed by a run of
    a DIFFERENT base). These matter because indel placement is inherently
    ambiguous inside them - the failure mode this replaces was
    hardcoding coordinates inferred from one sample's evidence.
  - Repeat-shift: right-anchoring an indel to its rightmost equivalent
    placement within a repeat/homopolymer run, the forensic 3'-shift
    convention.

Split out of the former repeat_regions.py (2026-09-23, user: "repeat_
regions.py contains a lot of functions that are not related to repeat
regions") - lh_bounds_pct/lh_bounds moved to their own heteroplasmy_
thresholds.py, a genuinely different concern (a reporting threshold, not
reference structure).
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
    report_separate_frame in process_fdstools_output_improved_better.py):
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


def load_reference_list(fasta_path):
    """Same sequence as load_reference above, as a mutable list of
    single-character bases instead of a string - needed wherever a
    caller bakes a substitution in place (e.g. process_mutect2_output_
    improved.py's major-SNP baking pre-pass does `reference[i] = v`,
    which a plain string can't support). Previously duplicated,
    byte-for-byte, as each of process_mutect2_output_improved.py's and
    process_fdstools_output_improved_better.py's own `load_reference`
    (one via Bio.SeqIO, one via manual line-parsing - verified to
    produce an identical sequence for the real reference file before
    consolidating here) (2026-09-23, user: "pull shared primitives").
    """
    return list(load_reference(fasta_path))


# --- Repeat-shift primitives -------------------------------------------
# Shared by process_mutect2_output_improved.py and process_fdstools_
# output_improved_better.py to right-anchor an indel to its rightmost
# equivalent placement within a repeat/homopolymer run - both VCF-row-
# oriented callers agree on the exact same algorithm (verified identical
# before consolidating), so one copy lives here instead of two.
#
# call_repeat_regions.py has its OWN, deliberately separate shift_
# insertion_right/shift_deletion_right - a simpler single-pass rotation
# used for direct BAM read/local-haplotype composition, a different
# context with different inputs (an observed read sequence, not a VCF
# REF/ALT pair) where the two-phase algorithm below doesn't apply the
# same way. Not unified with these: keeping a correctness-critical,
# unverified-equivalent algorithm change out of the (still real,
# standalone) BAM-direct tool was judged safer than a blind merge.

def rightmost_repeat_position(reference, pos, segment):
    pos -= 1
    while (
        pos + len(segment) < len(reference) and
        "".join(reference[pos + 1: pos + 1 + len(segment)]) == segment
    ):
        pos += len(segment)
    return pos + 1


def shift_insertion_right(reference, pos, segment):
    pos = rightmost_repeat_position(reference, pos, segment)
    for _ in range(len(segment) - 1):
        next_seq = reference[pos: pos + 1]
        if "".join(next_seq) == segment[0]:
            segment = segment[1:] + segment[0]
            pos += 1
        else:
            break
    return pos, segment


def shift_deletion_right(reference, pos, segment):
    pos = rightmost_repeat_position(reference, pos, segment) - len(segment) + 1
    for _ in range(len(segment) - 1):
        next_seq = reference[pos + len(segment) - 1: pos + len(segment)]
        if "".join(next_seq) == segment[0]:
            segment = segment[1:] + segment[0]
            pos += 1
        else:
            break
    return pos, segment


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
