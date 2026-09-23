#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import argparse
import sys
import os
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repeat_regions import find_homopolymer_runs, merge_runs, find_boundary_run_regions

# Utility: extract numeric position from sequence
def extract_position(seq):
    match = re.search(r"(\d+\.?\d*)", seq)
    return float(match.group(1)) if match else float('inf')

# Adjust for circular mtDNA positions (16570–16587 → 1–18)
def adjust_circular_position(pos):
    if 16570 <= pos <= 16587:
        return pos - 16569
    return pos

def load_reference(fasta_path):
    with open(fasta_path) as f:
        seq = "".join(line.strip() for line in f if not line.startswith(">"))
    return list(seq)

def rightmost_repeat_position(reference, pos, segment):
    pos -= 1
    while (
        pos + len(segment) < len(reference) and
        "".join(reference[pos + 1: pos + 1 + len(segment)]) == segment
    ):
        pos += len(segment)
    return pos + 1

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

LH_DROP_SENTINEL = "__BELOW_LH_FLOOR__"

def lh_bounds_pct(threshold_pct):
    """Symmetric floor/ceiling (percentage scale, 0-100) around a single
    length-heteroplasmy threshold, e.g. threshold_pct=10 -> (10, 90).
    Accepts either side (10 or 90) and always returns (floor, ceiling)
    with floor <= ceiling."""
    return min(threshold_pct, 100 - threshold_pct), max(threshold_pct, 100 - threshold_pct)

# IUPAC resolution for heteroplasmies
def resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold, IUPAC_CODES):
    seq = row['sequence']

    # A plain single-base deletion token gets the same lh_floor/lh_ceiling
    # major-minor branching the insertion-token case right below already
    # has - unconditionally replacing DEL with "-" regardless of
    # frequency was a real, pre-existing, genome-wide bug (not scoped to
    # anything this session touched): a 10.08%-frequency "A11038DEL" was
    # rendered as major "A11038-" instead of minor "A11038a", reading as
    # a spurious disagreement against MUTECT2's own correctly-formatted
    # "A11038a" at a near-identical frequency (2026-09-22, user, HG03366:
    # "in shared_frame what happened at A11038?"). Any other DEL shape
    # (a combined multi-position label, if one somehow still reaches this
    # point) falls back to the original unconditional behavior unchanged.
    if 'DEL' in seq:
        m = re.match(r'^([ACGT])(\d+)DEL$', seq)
        if m:
            lh_floor, lh_ceiling = lh_bounds_pct(length_heteroplasmy_threshold)
            if row['variant_frequency'] < lh_floor:
                return LH_DROP_SENTINEL
            ref, pos = m.group(1), m.group(2)
            if row['variant_frequency'] < lh_ceiling:
                return f"{ref}{pos}{ref.lower()}"
            return f"{ref}{pos}-"
        return seq.replace('DEL', '-')
    # Only a genuinely raw, unprocessed FDSTOOLS insertion token
    # ("16193.1C") gets the leading "-" and case treatment below - a
    # label some other step already fully formatted (e.g.
    # report_boundary_run's own "-16193.1c") also contains a "." but
    # must pass through unchanged, or it would get a second "-" prepended.
    if re.match(r'^\d+\.\d+[ACGT]$', seq):
        lh_floor, lh_ceiling = lh_bounds_pct(length_heteroplasmy_threshold)
        if row['variant_frequency'] < lh_floor:
            return LH_DROP_SENTINEL
        if row['variant_frequency'] < lh_ceiling:
            return '-' + seq[:-1] + seq[-1].lower()
        else:
            return '-' + seq
    if row['variant_frequency'] < 100 - min_variant_frequency_pct:
        match = re.match(r'([ACGT])(\d+)([ACGT])', seq)
        if match:
            ref, pos, alt = match.groups()
            code = IUPAC_CODES.get(frozenset([ref, alt]))
            if code:
                return f"{ref}{pos}{code}"
    return seq

# Boundary-run regions: a leading homopolymer run (e.g. the 16180-16183
# A-run) immediately followed by an extension run (e.g. the 16184-16193
# C-run) whose reference identity differs. Computed genome-wide from the
# reference alone (find_boundary_run_regions) - unlike the Mutect2/BAM
# side's equivalent fix, which stays deliberately scoped to just
# 16180-16193 (see BAM_OVERRIDE_REGIONS in process_mutect2_output_
# improved.py: Mutect2's own VCF rows are individual records, not
# complete per-read haplotypes, so there's no reliable per-region
# fallback signal to apply this genome-wide the way FDSTOOLS' complete
# per-read haplotype strings allow here). 16180-16193 itself is
# validated against s26-02989 and HG01799 (see report_boundary_run's
# docstring); the other 8 reference-derived regions this now covers
# follow the identical, already-validated shape and logic, but haven't
# each been checked against a real sample carrying a variant there.


def report_boundary_run(final, df_haplotypes, marker_total_reads, reference,
                         leading_start, leading_end, extension_start, extension_end,
                         lh_floor, lh_ceiling):
    """Report a leading homopolymer run and its immediately-adjacent
    extension run as two independent, position-relative counts, instead
    of the general reconciliation's single combined-and-frame-baked run.

    The general reconciliation above (see its own docstring) bakes every
    called substitution into one shared reference copy, then borrows
    THAT one frame to decide how every molecule's indels get spelled.
    That's only correct when the sample doesn't have real, differing
    subpopulations at the run's own boundary substitution - see the
    documented known limit in call_repeat_regions.py. s26-02989 is
    exactly that case: ~76% of reads carry A16182C, ~24% don't, and no
    single borrowed frame spells both groups correctly.

    This sidesteps the problem by never borrowing a frame at all:

    - The LEADING run (e.g. 16180-16183, the A-run) is reported as
      completely independent single positions - how many reads still
      show the reference base there, full stop. No repeat-shifting, no
      event-combining across positions. A position where every read
      still matches reference isn't reported at all; one with some but
      not all reads matching reference is reported as a minor call using
      the reference base's own letter (e.g. "A16182a"); one where
      essentially no read matches reference is reported as a major
      "REF+POS-" call (e.g. "A16183-").

    - The EXTENSION run (e.g. 16184-16193, the C-run) is reported as one
      shared cumulative "at least k extra bases" count - but a leading-
      run base converted AT THE BOUNDARY (e.g. A16182C, A16183C) counts
      as an extension exactly the same way a genuine insertion past the
      run's own far end does, since both are structurally
      indistinguishable once the run's own total length is what's being
      asked about. A deletion strictly INSIDE the extension run (e.g.
      T16189DEL) reduces that same total. Every read's own net extension
      is computed independently, straight from its own original
      haplotype string - no reference frame is baked or borrowed
      anywhere in this function.

    Validated read-for-read against s26-02989's real data (2026-09-16,
    106-read coverage): leading run -> A16180/A16181 no variant (106/106
    reference), A16182 "A16182a" 67.92% (72/106, i.e. NOT reference - same
    "lowercase = absent" convention as the general reconciliation's own
    deletion-direction rows, e.g. "C16193c"), A16183 "A16183-" 100.0%
    (0/106 reference); extension run -> "-16193.1c" 85.85% (91/106),
    "-16193.2c" 50.94% (54/106), "-16193.3c" 11.32% (12/106).
    """
    rows = []
    sub_re = re.compile(r'^([ACGT])(\d+)([ACGT])$')
    del_re = re.compile(r'^[ACGT](\d+)DEL$')
    ins_re = re.compile(r'^(\d+)\.(\d+)([ACGT])$')

    # ---- Leading run: independent per-position reference-presence ----
    # Counted from the RAW haplotype rows (df_haplotypes), not `final` -
    # final has already dropped whichever individual alt token didn't
    # clear the sample-wide min_vf threshold ON ITS OWN (e.g. a real but
    # rare G at a position whose C is the dominant alt), which would
    # silently undercount how often reference is actually absent here.
    # Every read genuinely showing something other than reference counts
    # toward loss regardless of how rare its own specific alt is - this
    # axis only asks "is reference here or not", not "which alt".
    leading_token_re = re.compile(r'^([ACGT])(\d+)(?:[ACGT]|DEL)$')
    non_ref_by_marker_pos = {}
    for _, row in df_haplotypes.iterrows():
        seq_str = str(row["sequence"]).strip()
        if seq_str == "Other sequences":
            continue
        marker, n = row["marker"], row["total"]
        seen = set()
        for tok in seq_str.split():
            m = leading_token_re.match(tok)
            if not m:
                continue
            p = int(m.group(2))
            if leading_start <= p <= leading_end and p not in seen:
                non_ref_by_marker_pos[(marker, p)] = non_ref_by_marker_pos.get((marker, p), 0) + n
                seen.add(p)

    for pos in range(leading_start, leading_end + 1):
        ref_here = reference[pos - 1].upper()
        for marker, coverage in marker_total_reads.items():
            non_ref = non_ref_by_marker_pos.get((marker, pos), 0)
            if not coverage or not non_ref:
                continue
            ref_total = coverage - non_ref
            loss_pct = round(100.0 * non_ref / coverage, 2)
            if loss_pct < lh_floor:
                continue
            # Both branches' own explanatory variant_note text is
            # suppressed - not deleted, kept here so it can be turned
            # back on - per the user (2026-09-22): "too much information
            # that is difficult to follow ... I don't want to see it" /
            # "let's not delete it" / "suppress it". The FMP label,
            # total, and frequency below are untouched either way.
            # major_note = (f"reference {ref_here} observed at only {int(ref_total)}/{coverage} reads here "
            #               f"- the rest are reported as part of the {extension_start}-{extension_end} run below")
            # minor_note = (f"reference {ref_here} retained at {int(ref_total)}/{coverage} reads here "
            #               f"- the rest are reported as part of the {extension_start}-{extension_end} run below")
            if loss_pct >= lh_ceiling:
                rows.append({
                    "sequence": f"{ref_here}{pos}-",
                    "total": coverage - ref_total,
                    "interpolated_total_coverage": coverage,
                    "is_noise_or_low_frq": False,
                    "num_markers": 1,
                    "variant_frequency": loss_pct,
                    "marker": marker,
                    "position": float(pos),
                })
            else:
                # Same quantity as the major case above (loss_pct - the
                # fraction WITHOUT reference here) for both total and
                # variant_frequency, matching the existing deletion-
                # direction convention elsewhere in this file (e.g.
                # "C16193c"): the lowercase letter always reports how
                # often that base is ABSENT, not how often it's observed
                # - only the major/minor case of the suffix differs.
                rows.append({
                    "sequence": f"{ref_here}{pos}{ref_here.lower()}",
                    "total": coverage - ref_total,
                    "interpolated_total_coverage": coverage,
                    "is_noise_or_low_frq": False,
                    "num_markers": 1,
                    "variant_frequency": loss_pct,
                    "marker": marker,
                    "position": float(pos),
                })

    # ---- Extension run: unified left+right cumulative "at least k" ----
    extension_base = reference[extension_start - 1].upper()
    net_by_marker = {}
    for _, row in df_haplotypes.iterrows():
        seq_str = str(row["sequence"]).strip()
        if seq_str == "Other sequences":
            continue
        tokens = seq_str.split()
        marker, n = row["marker"], row["total"]

        subbed = {int(m.group(2)): m.group(3) for tok in tokens for m in [sub_re.match(tok)] if m}
        left = 0
        p = extension_start - 1
        while p >= leading_start and subbed.get(p) == extension_base:
            left += 1
            p -= 1

        # Mirror of `left`, walking the other direction: a substitution
        # just past the run's own far end (e.g. "A16194C", reference A at
        # 16194 substituted to C) extends the run rightward on that read
        # exactly the way a leading-run substitution extends it leftward -
        # found in real data (HG02389: "A16194C 16194.1C" together, a
        # substitution immediately followed by an insertion anchored at
        # the position IT creates, not at the reference's own extension_end).
        # Self-limiting (stops at the first non-matching position), no
        # separate upper bound needed - there's no far-side region boundary
        # to stop at the way leading_start bounds the left-hand walk.
        right_sub = 0
        p = extension_end + 1
        while subbed.get(p) == extension_base:
            right_sub += 1
            p += 1
        right_edge = extension_end + right_sub

        deleted = {int(m.group(1)) for tok in tokens for m in [del_re.match(tok)] if m}
        internal = sum(1 for d in deleted if extension_start <= d <= extension_end)

        # Mirror of `right`'s own right_edge anchoring, not "anywhere in
        # the leading run": an insertion only genuinely extends the run
        # when it sits immediately adjacent to wherever the run's own
        # (possibly substitution-extended) left edge currently is for
        # THIS read - extension_start - 1 - left, exactly mirroring
        # right_edge. An insertion anchored any earlier still has an
        # untouched reference base of the leading run's own identity
        # sitting between it and the run itself, so it can't be folded in
        # the same way (real case, HG01799: "16182.1C" sits between the
        # 3rd and 4th A of a plain AAAA leading run with left=0, so
        # left_edge is 16183, not 16182 - treating it as "one more C" was
        # wrong, since the 4th A is still right there between the
        # insertion and the actual C run; 2026-09-22, user: "the problem
        # here is that the insertion of 16182.1C occurs before the fourth
        # A and so we cannot shift it to 16193"). Such a stray insertion
        # isn't dropped from `final` by the caller either now, so it
        # passes through and reports itself, at its own true anchor.
        left_edge = extension_start - 1 - left
        leading_insertion = 0
        for tok in tokens:
            m = ins_re.match(tok)
            if m and m.group(3) == extension_base and int(m.group(1)) == left_edge:
                leading_insertion = max(leading_insertion, int(m.group(2)))

        # Anchored at right_edge, not a fixed extension_end - once a
        # right_sub walk has extended the run to e.g. 16194 for this read,
        # a further insertion continues from there ("16194.1C"), not from
        # the reference's own unextended boundary.
        right = 0
        for tok in tokens:
            m = ins_re.match(tok)
            if m and int(m.group(1)) == right_edge:
                right = max(right, int(m.group(2)))

        net = left + leading_insertion - internal + right_sub + right
        net_by_marker.setdefault(marker, []).append((net, n))

    for marker, entries in net_by_marker.items():
        coverage = marker_total_reads.get(marker)
        if not coverage:
            continue
        k = 1
        while True:
            reads = sum(n for net, n in entries if net >= k)
            if not reads:
                break
            pct = round(reads / coverage * 100, 2)
            if pct < lh_floor:
                break
            is_major = pct >= lh_ceiling
            base = extension_base if is_major else extension_base.lower()
            # Explanatory variant_note suppressed - not deleted, kept
            # here for the same reason as the leading-run notes above.
            # note = (f"{reads} of {coverage} reads ({pct}%) have at least {k} extra "
            #         f"{extension_base}{'s' if k != 1 else ''} in the {extension_start}-{extension_end} "
            #         f"run - from substitutions at the run's own leading edge "
            #         f"({leading_start}-{leading_end}), insertions past {extension_end}, or both combined")
            rows.append({
                "sequence": f"-{extension_end}.{k}{base}",
                "total": reads,
                "interpolated_total_coverage": coverage,
                "is_noise_or_low_frq": False,
                "num_markers": 1,
                "variant_frequency": pct,
                "marker": marker,
                "position": float(extension_end),
            })
            k += 1

    return rows


def dominant_boundary_molecule(df_haplotypes, marker_total_reads, reference,
                                leading_start, leading_end, extension_start, extension_end,
                                lh_floor, lh_ceiling):
    """ISFG/EMPOP-style reporting for a boundary-run region: report only the
    single most common molecule and only if it differs from reference, no
    lowercase/minor calls at all. This is what W. Parson's own reading of
    the DNA Commission guidelines calls for (Recommendation #11's dominant-
    type convention; his note, 2026-09-21, that the lowercase mixture
    convention was intended for point heteroplasmy and indel mixtures
    outside homopolymer runs, not partial states within one) - an
    alternative to report_boundary_run's own two-axis scheme above, not a
    replacement for it (see homopolymer_reporting mode in
    process_fdstools_sast).

    Groups df_haplotypes' own real per-read (leading pattern, interrupt-
    position confirm, extension net) combinations, per marker, and
    reports whichever exact combination has the most reads - not each
    axis's own separate marginal majority, which need not correspond to
    any single real read (see plot_boundary_run.py's collect_patterns,
    validated against HG01799 and s26-02989 2026-09-21/22: HG01799's
    dominant molecule matches reference exactly, 518/828 reads, 62.6% -
    nothing would be reported for it under this mode; s26-02989's
    dominant molecule loses both leading positions and carries 2 extra
    extension bases, 33/108 reads, only 30.6%, so most reads still don't
    match even the dominant call there).

    The extension run's own tolerated single-base interruption (e.g.
    16184-16193's T16189, a plain point substitution when confirmed as
    C - see find_boundary_run_regions) joins the same joint key, not left
    to leak through the separate, generic point-substitution pipeline as
    its own differently-shaped row: the whole point of dominant_only is
    one unified report per region, not the dominant axis's own read split
    across mechanisms (2026-09-22, user, after seeing exactly that split
    for this region: "we should also report only the dominant molecule in
    dominant_only_boundary_run for this region").
    """
    sub_re = re.compile(r'^([ACGT])(\d+)([ACGT])$')
    del_re = re.compile(r'^[ACGT](\d+)DEL$')
    ins_re = re.compile(r'^(\d+)\.(\d+)([ACGT])$')
    leading_token_re = re.compile(r'^([ACGT])(\d+)(?:[ACGT]|DEL)$')
    extension_base = reference[extension_start - 1].upper()
    interrupt_positions = [
        p for p in range(extension_start, extension_end + 1)
        if reference[p - 1].upper() != extension_base
    ]

    def read_pattern(seq_str):
        tokens = seq_str.split()
        non_ref = set()
        for tok in tokens:
            m = leading_token_re.match(tok)
            if m:
                p = int(m.group(2))
                if leading_start <= p <= leading_end:
                    non_ref.add(p)
        leading_presence = tuple(p not in non_ref for p in range(leading_start, leading_end + 1))

        subbed = {int(m.group(2)): m.group(3) for tok in tokens for m in [sub_re.match(tok)] if m}
        deleted = {int(m.group(1)) for tok in tokens for m in [del_re.match(tok)] if m}
        # An interrupt position reads as "confirmed" (reported as a plain
        # substitution, e.g. "T16189C") whether it's a literal
        # substitution token or an explicit deletion - a T sitting
        # directly between two stretches of the run's own base has no
        # other difference from reference around it, so deleting it and
        # substituting it to the run's base are sequence-identical
        # (closing the T-shaped gap either way merges the flanking C's
        # into one run); this matches the row's own existing FMP label,
        # which already reports it as "T16189C" via the pre-existing
        # fold-in convention regardless of which token type produced it.
        # Reporting it as "confirmed" without ALSO accounting for the
        # length axis would overstate the run by one base, though - a
        # literal deletion doesn't actually add a C the way a real
        # substitution does, so internal (below) still counts it, letting
        # net correctly come out negative and trigger a real compensating
        # C-deletion elsewhere in the run (2026-09-22, user, after seeing
        # "T16189-" sit right next to the row's own "T16189C" label: "if
        # we report T16189C then there should be a deletion of C as
        # well, no?").
        def interrupt_state(p):
            return "confirmed" if (subbed.get(p) == extension_base or p in deleted) else "reference"
        interrupt_states = tuple(interrupt_state(p) for p in interrupt_positions)

        left = 0
        p = extension_start - 1
        while p >= leading_start and subbed.get(p) == extension_base:
            left += 1
            p -= 1
        # Mirror of `left` - see report_boundary_run's own copy of this
        # logic above for the full rationale (HG02389's real "A16194C
        # 16194.1C").
        right_sub = 0
        p = extension_end + 1
        while subbed.get(p) == extension_base:
            right_sub += 1
            p += 1
        right_edge = extension_end + right_sub
        # Includes interrupt_positions' own deletions again - a literal
        # deletion there is reported as "confirmed" above (matching the
        # row's own "T16189C" label), but it doesn't actually add a C the
        # way a real substitution would, so it still has to count against
        # the run's true length here, or the run would be overstated by
        # one base with nothing to correct it.
        internal = sum(1 for d in deleted if extension_start <= d <= extension_end)
        # Mirror of right_edge, not "anywhere in the leading run" - see
        # report_boundary_run's own copy of this same fix above for the
        # full rationale (HG01799's real "16182.1C", sitting between the
        # 3rd and 4th A of a plain AAAA run, wrongly treated as extending
        # the C run when the 4th A is still right there in between).
        left_edge = extension_start - 1 - left
        leading_insertion = 0
        for tok in tokens:
            m = ins_re.match(tok)
            if m and m.group(3) == extension_base and int(m.group(1)) == left_edge:
                leading_insertion = max(leading_insertion, int(m.group(2)))
        # A leading-run insertion anchored anywhere OTHER than left_edge
        # (e.g. HG01799's real "16182.1C", with the 4th A still sitting
        # between it and the C run) doesn't extend the run - but it's
        # still a real difference from reference and belongs in the joint
        # key just like interrupt_states does, or it silently vanishes
        # from the dominant molecule's own report even when it's present
        # in virtually every read (2026-09-22, user: "-16182.1C also
        # belongs to the dominant molecule so why is it not in variant
        # note? every variant between 16180 and 16193 (plus insertions)
        # should be in variant note if it's part of a dominant molecule").
        stray_by_anchor = {}
        for tok in tokens:
            m = ins_re.match(tok)
            if (m and m.group(3) == extension_base
                    and leading_start <= int(m.group(1)) <= leading_end
                    and int(m.group(1)) != left_edge):
                a = int(m.group(1))
                stray_by_anchor[a] = max(stray_by_anchor.get(a, 0), int(m.group(2)))
        stray_insertions = tuple(sorted(stray_by_anchor.items()))
        right = 0
        for tok in tokens:
            m = ins_re.match(tok)
            if m and int(m.group(1)) == right_edge:
                right = max(right, int(m.group(2)))
        net = left + leading_insertion - internal + right_sub + right
        return leading_presence, interrupt_states, stray_insertions, net

    counts_by_marker = {}
    for _, row in df_haplotypes.iterrows():
        seq_str = str(row["sequence"]).strip()
        if seq_str == "Other sequences":
            continue
        marker, n = row["marker"], row["total"]
        pattern = read_pattern(seq_str)
        d = counts_by_marker.setdefault(marker, {})
        d[pattern] = d.get(pattern, 0) + n

    rows = []
    for marker, coverage in marker_total_reads.items():
        patterns = counts_by_marker.get(marker)
        if not patterns or not coverage:
            continue
        (leading_presence, interrupt_states, stray_insertions, net), count = max(patterns.items(), key=lambda kv: kv[1])
        pct = round(100.0 * count / coverage, 2)
        if pct < lh_floor:
            continue
        note = f"dominant molecule ({int(count)}/{int(coverage)} reads, {pct}%)"

        if (all(leading_presence) and net == 0 and not stray_insertions
                and all(s == "reference" for s in interrupt_states)):
            continue  # dominant molecule matches reference - nothing to report

        for anchor, ins_count in stray_insertions:
            for k in range(1, ins_count + 1):
                rows.append({
                    "sequence": f"-{anchor}.{k}{extension_base}", "total": count,
                    "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                    "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                    "marker": marker, "position": float(anchor),
                })
        for offset, present in enumerate(leading_presence):
            if present:
                continue
            pos = leading_start + offset
            ref_here = reference[pos - 1].upper()
            rows.append({
                "sequence": f"{ref_here}{pos}-", "total": count,
                "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                "marker": marker, "position": float(pos),
            })
        for p, state in zip(interrupt_positions, interrupt_states):
            if state == "reference":
                continue
            ref_here = reference[p - 1].upper()
            rows.append({
                "sequence": f"{ref_here}{p}{extension_base}", "total": count,
                "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                "marker": marker, "position": float(p),
            })
        if net > 0:
            for k in range(1, net + 1):
                rows.append({
                    "sequence": f"-{extension_end}.{k}{extension_base}", "total": count,
                    "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                    "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                    "marker": marker, "position": float(extension_end),
                })
        elif net < 0:
            # A compensating loss elsewhere in the run - most often from
            # an interrupt position reported as "confirmed" above (e.g.
            # "T16189C") without a real substitution behind it, so the
            # nominal +1 that implies has to be given back somewhere (
            # 2026-09-22, user: "if we report T16189C then there should
            # be a deletion of C as well, no?"), but not exclusively - any
            # other uncompensated deletion within the extension run lands
            # here too. Same shift-right convention as dominant_shared_
            # frame_molecule's own net<0 branch and the general
            # reconciliation's deletion-shift machinery elsewhere in this
            # file, so it reads the same way shared_frame already does
            # for this exact molecule (e.g. "C16193-").
            segment = extension_base * abs(net)
            shifted_pos, shifted_segment = shift_deletion_right(reference, extension_end - abs(net) + 1, segment)
            for i, base in enumerate(shifted_segment):
                pos = shifted_pos + i
                rows.append({
                    "sequence": f"{base}{pos}-", "total": count,
                    "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                    "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                    "marker": marker, "position": float(pos),
                })
    return rows


def find_baked_positions_by_run(reference, homopolymer_runs):
    """Every position within a homopolymer run where reference_baked
    differs from the true, unbaked reference (see dominant_shared_frame_
    molecule below) - purely a property of the reference and the run's
    own boundaries, not of any sample's actual reads.

    Baking has no majority threshold, so this includes both a genuinely
    fixed substitution (T16189C at 100%, s26-02989) and a real point-
    heteroplasmy (A16182C at only 60%, same sample - 35 of 108 reads
    genuinely keep the true reference A there) alike. Both fold into the
    same joint dominant-molecule report below now, not just the
    heteroplasmic one: reporting only A16182C's own confirm state left
    T16189C and A16183C to leak through the separate, generic point-
    substitution/IUPAC pipeline as their own, differently-shaped rows,
    which is exactly the fragmentation the user asked to close (2026-
    09-22): "we want this for the whole region from 16180 to 16193
    (including insertions afterwards). I think FDSTOOLS already does
    this, we just need to identify the most common molecule from sast
    file and report it." A position's own raw point-substitution/
    deletion row gets dropped from `final` by the caller for the same
    reason - one report per position, not two.
    """
    baked_by_run = {}
    for run_start, run_end, run_base, _length in homopolymer_runs:
        baked = [p for p in range(run_start, run_end + 1) if reference[p - 1].upper() != run_base]
        if baked:
            baked_by_run[(run_start, run_end, run_base)] = baked
    return baked_by_run


def dominant_shared_frame_molecule(df_haplotypes, marker_total_reads, reference, reference_baked,
                                    homopolymer_runs, baked_positions_by_run, lh_floor, lh_ceiling):
    """ISFG/EMPOP-style dominant-molecule reporting (see dominant_boundary_
    molecule's own docstring for the full rationale and citation), built
    on the shared_frame model instead of boundary_run's two-axis one -
    the only option for a plain homopolymer run with no different-base
    leading run to split off (e.g. chrM 303-315, the "310" region), and a
    genuine second, independent determination for boundary-run-shaped
    ones too, worth comparing against dominant_only_boundary_run's own
    answer rather than assumed to always agree with it. Generic over
    whatever `homopolymer_runs` it's handed - process_fdstools_sast
    restricts that list to target_regions (see there), which now
    includes 303-315 alongside chrM 16180's own region.

    Per (marker, individual homopolymer run - not the merged, possibly
    mixed-base regions used for the deletion-shift machinery above, an
    insertion or deletion's own base has to match the ONE run it affects),
    each read's net length is computed directly from its own raw
    haplotype tokens in one pass - insertions of the run's own base
    (inside it, or past either edge via a substitution extending into
    that edge - same per-read math as the shared_frame insertion loop
    above) minus deletions within the run - and whichever net value has
    the most reads is reported, only if it differs from reference.
    Deliberately NOT built from the separately-extracted `events`/
    `ins_events_by_run` structures above: those are pre-aggregated across
    reads for the full-distribution report and would reintroduce exactly
    the marginal-vs-joint mistake dominant_boundary_molecule's own
    docstring warns against - a read with both an insertion and a
    deletion needs its own single net computed together, not two
    separately-aggregated axes recombined after the fact.

    Alongside net length, each read's own confirm/revert state at every
    one of that run's baked_positions_by_run entries (from
    find_baked_positions_by_run above - every baked position, fixed or
    heteroplasmic alike, not just the heteroplasmic ones) joins the same
    grouping key - the read's true joint molecule across the whole
    region, not net and each point decided independently. The winning
    group's confirmed positions are reported as plain substitutions (the
    caller drops their now-superseded raw rows from `final` for this same
    reason: one unified report per region, not one row per position
    scattered across different mechanisms).
    """
    sub_re = re.compile(r'^([ACGT])(\d+)([ACGT])$')
    del_re = re.compile(r'^[ACGT](\d+)DEL$')
    ins_re = re.compile(r'^(\d+)\.(\d+)([ACGT])$')

    counts_by_run = {}
    for _, row in df_haplotypes.iterrows():
        seq_str = str(row["sequence"]).strip()
        if seq_str == "Other sequences":
            continue
        tokens = seq_str.split()
        marker, n = row["marker"], row["total"]
        subbed = {int(m.group(2)): m.group(3) for tok in tokens for m in [sub_re.match(tok)] if m}
        deleted = {int(m.group(1)) for tok in tokens for m in [del_re.match(tok)] if m}

        for run_start, run_end, run_base, _length in homopolymer_runs:
            ins_by_anchor = {}
            for tok in tokens:
                m = ins_re.match(tok)
                if m and m.group(3) == run_base and run_start <= int(m.group(1)) <= run_end:
                    a = int(m.group(1))
                    ins_by_anchor[a] = max(ins_by_anchor.get(a, 0), int(m.group(2)))
            inside = sum(ins_by_anchor.values())

            left_sub = 0
            p = run_start - 1
            while subbed.get(p) == run_base:
                left_sub += 1
                p -= 1
            left_edge = run_start - left_sub

            right_sub = 0
            p = run_end + 1
            while subbed.get(p) == run_base:
                right_sub += 1
                p += 1
            right_edge = run_end + right_sub

            left_edge_ins = 0
            if left_sub > 0:
                for tok in tokens:
                    m = ins_re.match(tok)
                    if m and int(m.group(1)) == left_edge - 1:
                        left_edge_ins = max(left_edge_ins, int(m.group(2)))
            right_edge_ins = 0
            if right_sub > 0:
                for tok in tokens:
                    m = ins_re.match(tok)
                    if m and int(m.group(1)) == right_edge:
                        right_edge_ins = max(right_edge_ins, int(m.group(2)))

            # An insertion of the run's own base anchored just short of
            # either edge (still a real reference base of some OTHER
            # identity sitting between it and the run) doesn't extend the
            # run either - same mirror fix as dominant_boundary_molecule's
            # own left_edge/stray_insertions above (2026-09-22, user,
            # after HG01799's real "16182.1C": "also for shared_frame").
            # Windowed to 4 positions either side, matching the min_length
            # =4 homopolymer threshold used throughout this file - this
            # function has no separate leading_start/leading_end of its
            # own to bound the search with, unlike dominant_boundary_
            # molecule's region-specific one.
            stray_by_anchor = {}
            left_excluded = left_edge - 1 if left_sub > 0 else None
            right_excluded = right_edge if right_sub > 0 else None
            for tok in tokens:
                m = ins_re.match(tok)
                if not (m and m.group(3) == run_base):
                    continue
                a = int(m.group(1))
                in_left_window = run_start - 4 <= a < run_start and a != left_excluded
                in_right_window = run_end < a <= run_end + 4 and a != right_excluded
                if in_left_window or in_right_window:
                    stray_by_anchor[a] = max(stray_by_anchor.get(a, 0), int(m.group(2)))
            stray_insertions = tuple(sorted(stray_by_anchor.items()))

            deletion_total = sum(1 for d in deleted if run_start <= d <= run_end)

            net = inside + left_sub + right_sub + left_edge_ins + right_edge_ins - deletion_total
            key = (marker, run_start, run_end, run_base)
            baked = baked_positions_by_run.get((run_start, run_end, run_base), [])
            # "Confirmed" (reported as a plain substitution) on a literal
            # deletion too, not just a literal substitution token - same
            # fix as dominant_boundary_molecule's own interrupt_state
            # above, same rationale: a baked position deleted with
            # nothing else different around it is sequence-identical to
            # substituting it to the run's base and losing one base
            # elsewhere in the run, and deletion_total above already
            # counts it either way, so the compensating loss still shows
            # up correctly via net (2026-09-22, user: "why in the
            # shared_frame the T16189C is not reported as part of the
            # dominant molecule? in the variant note").
            confirm = tuple(subbed.get(p) == run_base or p in deleted for p in baked)
            d = counts_by_run.setdefault(key, {})
            d[(confirm, stray_insertions, net)] = d.get((confirm, stray_insertions, net), 0) + n

    rows = []
    for (marker, run_start, run_end, run_base), pattern_counts in counts_by_run.items():
        coverage = marker_total_reads.get(marker, 0)
        if not coverage:
            continue
        (confirm, stray_insertions, net), count = max(pattern_counts.items(), key=lambda kv: kv[1])
        baked = baked_positions_by_run.get((run_start, run_end, run_base), [])
        if net == 0 and not any(confirm) and not stray_insertions:
            continue  # dominant molecule matches reference - nothing to report
        pct = round(100.0 * count / coverage, 2)
        if pct < lh_floor:
            continue
        note = f"dominant molecule ({int(count)}/{int(coverage)} reads, {pct}%)"

        for anchor, ins_count in stray_insertions:
            for k in range(1, ins_count + 1):
                rows.append({
                    "sequence": f"-{anchor}.{k}{run_base}", "total": count,
                    "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                    "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                    "marker": marker, "position": float(anchor),
                })

        for p, is_confirmed in zip(baked, confirm):
            if is_confirmed:
                rows.append({
                    "sequence": f"{reference[p - 1].upper()}{p}{run_base}", "total": count,
                    "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                    "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                    "marker": marker, "position": float(p),
                })

        if net > 0:
            for k in range(1, net + 1):
                rows.append({
                    "sequence": f"-{run_end}.{k}{run_base}", "total": count,
                    "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                    "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                    "marker": marker, "position": float(run_end),
                })
        elif net < 0:
            segment = run_base * abs(net)
            shifted_pos, shifted_segment = shift_deletion_right(reference_baked, run_end - abs(net) + 1, segment)
            for i, base in enumerate(shifted_segment):
                pos = shifted_pos + i
                rows.append({
                    "sequence": f"{base}{pos}-", "total": count,
                    "interpolated_total_coverage": coverage, "is_noise_or_low_frq": False,
                    "num_markers": 1, "variant_frequency": pct, "variant_note": note,
                    "marker": marker, "position": float(pos),
                })
    return rows


def load_marker_ranges(filepath):
    marker_ranges = {}
    in_position_block = False
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("[genome_position]"):
                in_position_block = True
                continue
            if line.startswith("[") and in_position_block:
                break
            if in_position_block and "=" in line:
                marker, values = line.split("=")
                parts = [v.strip() for v in values.split(",")]
                if len(parts) >= 3:
                    chrom, start, end = parts[:3]
                    marker_ranges[marker.strip()] = f"{start}-{end}"
    return marker_ranges

# Main processing function
def process_fdstools_sast(file_path, marker_map_path, output_file, reference_fasta, min_variant_frequency_pct=5.0, depth_threshold=10, length_heteroplasmy_threshold=90.0, homopolymer_reporting="boundary_run"):
    reference = load_reference(reference_fasta)
    IUPAC_CODES = {
        frozenset(["A", "G"]): "R", frozenset(["C", "T"]): "Y",
        frozenset(["A", "C"]): "M", frozenset(["G", "T"]): "K",
        frozenset(["G", "C"]): "S", frozenset(["A", "T"]): "W"
    }

    marker_map = load_marker_ranges(marker_map_path)

    df = pd.read_csv(file_path, sep="\t", dtype=str)
    df = df.drop(columns=[
        "total_mp_max", "forward_pct", "forward", "forward_mp_sum",
        "forward_mp_max", "reverse", "reverse_mp_sum", "reverse_mp_max"
    ], errors="ignore")

    df["total_mp_sum"] = pd.to_numeric(df["total_mp_sum"], errors="coerce").fillna(0)
    df["total"] = pd.to_numeric(df["total"], errors="coerce").fillna(0)

    # Flag every known amplicon (from marker_map, not just what happens
    # to appear in this file) whose TOTAL observed reads - summed across
    # every sequence FDSTOOLS reported for it, including its own
    # low-frequency "Other sequences" bucket - fall below the depth
    # threshold. Checking only markers reported as a single row (the
    # previous version) misses two real failure modes: a low-depth
    # amplicon reported as several small fragments (summed total still
    # below threshold, but marker_counts > 1 so it was never even
    # looked at, and each fragment then just quietly gets filtered out
    # downstream with no visible trace), and a complete dropout (zero
    # reads, so the marker never appears as a row at all - can't be
    # caught by inspecting df's own rows, only by checking against the
    # full known marker list).
    marker_total_reads = df.groupby("marker")["total"].sum()
    single_low_coverage = pd.DataFrame([
        {"marker": marker, "total": marker_total_reads.get(marker, 0), "sequence": "LOW"}
        for marker in marker_map
        if marker_total_reads.get(marker, 0) < depth_threshold
    ])


    df["total"] = df["total"].fillna(0)
    df["total_mp_sum"] = df["total_mp_sum"].fillna(0)

    # Snapshot before the explode below reassigns df to one row per
    # token - the repeat-region reconciliation further down needs the
    # original, space-joined per-haplotype sequence strings intact.
    df_haplotypes = df.copy()

    # Exact per-marker coverage (excludes "Other sequences", no
    # percentage math at all) - used directly as interpolated_total_
    # coverage below, replacing the old ceil(total/total_mp_sum%) back-
    # calculation, which reconstructed coverage from FDSTOOLS' own
    # total_mp_sum field rounded to just 1 decimal place and could
    # overshoot the true depth by 1 whenever the true percentage sat
    # just past a rounding boundary (2026-09-22, user, s24-12883a: a
    # real 19-read haplotype whose FDSTOOLS-reported 82.6% - truncated
    # from the true 82.6087% - back-calculated to a coverage of 20, not
    # 19). Computed from df_haplotypes (before the explode below splits
    # "Other sequences" into separate "Other"/"sequences" tokens) - the
    # same computation the homopolymer-region reconciliation further
    # down already relies on for exactly this reason.
    is_other = df_haplotypes["sequence"].astype(str).str.strip() == "Other sequences"
    marker_total_reads = df_haplotypes[~is_other].groupby("marker")["total"].sum()

    # Step 5: Split multiple variants
    df = df.assign(sequence=df["sequence"].str.split())
    df = df.explode("sequence").reset_index(drop=True)

    df["interpolated_total_coverage"] = df["marker"].map(marker_total_reads).fillna(0).astype("Int64")

    grouped = df.groupby(["marker", "sequence"], as_index=False).agg(
        total=("total", "sum"),
        interpolated_total_coverage=("interpolated_total_coverage", "max"),
    )
    grouped["interpolated_total_coverage"] = grouped["interpolated_total_coverage"].clip(lower=1)

    final = grouped.groupby("sequence", as_index=False).agg(
        marker=("marker", "first"),
        total=("total", "sum"),
        interpolated_total_coverage=("interpolated_total_coverage", "sum"),
        num_markers=("marker", "nunique")
    )

    final["variant_frequency"] = (final["total"] / final["interpolated_total_coverage"] * 100).round(2)
    # final["variant_frequency_wo_noise_or_low_frq"] = (final["total"] / final["total_wo_noise_or_low_frq"] * 100).round(2)
    final["position"] = final["sequence"].apply(extract_position)

    drop_seqs = ["Other", "sequences", "REF", "N3107DEL", "No", "data"]
    final = final[(~final["sequence"].isin(drop_seqs))]

    final["is_noise_or_low_frq"] = (final["sequence"].isin(["Other sequences"])) | (final["variant_frequency"] < min_variant_frequency_pct)
    final = final[~final["is_noise_or_low_frq"]]

    # --- Repeat-region deletion/substitution reconciliation ---
    #
    # FDSTOOLS already resolves each read into one definite haplotype
    # string (no alignment-placement ambiguity the way a raw aligner
    # CIGAR has inside a homopolymer), so every token in an original
    # haplotype row is a trustworthy observation. Two things still need
    # reconciling before the exploded per-position rows are report-ready:
    #
    # 1) A multi-base deletion explodes into several single-position DEL
    #    rows. Which of those actually belong to one atomic event can
    #    only be read off the ORIGINAL, pre-explosion haplotype string -
    #    e.g. "T310DEL" (11 reads) and "T310DEL C311DEL" (13 reads) are
    #    two genuinely different, disjoint read populations (1bp vs 2bp
    #    short) that happen to nest at adjacent positions, whereas
    #    "A16182DEL A16183DEL" (5 reads) is one atomic 2bp event. Trying
    #    to reconstruct this after explosion+re-aggregation (matching
    #    totals across positions) is ambiguous; going back to the source
    #    haplotype rows is not.
    #
    # 2) A read with a deletion at a position has no base call there at
    #    all - it can't be "the substitution" or "not the substitution",
    #    it's simply uninformative for that axis. So a co-located
    #    substitution's frequency must exclude those reads from its
    #    denominator entirely (not just its numerator), regardless of
    #    whether the deletion itself individually clears the length-
    #    heteroplasmy floor: a below-floor deletion isn't trusted as its
    #    own reportable event, but the reads still plainly have no base
    #    call at that position, and leaving them in the denominator
    #    makes a homoplasmic substitution look artificially heteroplasmic.
    #
    # The repeat-shift itself (shift_deletion_right) is only mathematically
    # valid for a single-base homopolymer run (it slides by checking
    # "does the next reference base match the first base of what's being
    # deleted") - so this whole reconciliation is scoped to deletion
    # events fully contained in a reference-derived homopolymer region.
    # A multi-base tandem repeat (e.g. the 8281-8289 9bp CA-repeat
    # deletion, bases C,C,C,C,C,T,C,T,A) isn't one - applying this logic
    # there produced inconsistent, sometimes-overlapping shifted anchors.
    # Anything outside a homopolymer region falls back to the simpler,
    # single-pair naive-position merge this replaces for that case.

    sub_pattern = re.compile(r'^([ACGT])(\d+)([ACGT])$')
    del_token_re = re.compile(r'^[ACGT](\d+)DEL$')

    # Bake every called substitution into one shared reference copy, so
    # an indel's rightward shift search sees the true, fully-substituted
    # run - not just whichever single substitution sits at its own naive
    # position. (A per-event-local bake would stop the search short:
    # e.g. a 16183 deletion's search would stop at 16188 without also
    # knowing 16189 is substituted to C.)
    reference_baked = list(reference)
    for _, row in final.iterrows():
        m = sub_pattern.match(str(row["sequence"]))
        if m:
            reference_baked[int(m.group(2)) - 1] = m.group(3)

    # Every piece of homopolymer-region reporting below - the deletion
    # repeat-shift machinery that predates homopolymer_reporting as a
    # parameter, and this session's own boundary-run dispatch and general-
    # path insertion reconciliation - is restricted to just two regions,
    # collected into target_regions: chrM 16180 through wherever the
    # poly-C tract it feeds into actually ends (16189 T>C's own famous
    # consequence), and chrM 303-315 (7 C's, T, 5 C's - the same
    # T-interrupted-run shape, but with no different-base leading run
    # long enough to split off: chrM 300-302 is only 3 A's, below the
    # min_length=4 homopolymer threshold used throughout this file, so
    # find_boundary_run_regions doesn't find this one at all - identified
    # via the same merge_runs(gap=1) convention homopolymer_regions
    # itself already uses below instead. 2026-09-23, user, s25-04643's
    # real T>C at 310: "can we try extending our fdstools shared frame
    # and boundary run to this region as well? ... the three As are not
    # considered a homopolymer stretch, right?"). Every other boundary-
    # run region find_boundary_run_regions detects genome-wide (353, 452,
    # 1527, 5488, 5597, 6416, 14531, 14616), plus every other plain 4+
    # run, is left completely untouched - passed through exactly as
    # FDSTOOLS reported it - for every homopolymer_reporting mode.
    _target_br = next(
        (br for br in find_boundary_run_regions("".join(reference).upper()) if br["leading"][0] == 16180),
        None
    )
    _seed_regions = [(_target_br["leading"][0], _target_br["extension"][1])] if _target_br else []
    _merged_303 = next(
        (r for r in merge_runs(find_homopolymer_runs("".join(reference).upper(), 4), 1) if r[0] <= 303 <= r[1]),
        None
    )
    if _merged_303:
        _seed_regions.append((_merged_303[0], _merged_303[1]))

    # Each seed gets independently widened with whatever reference_baked's
    # own run actually reaches, in case baking (no majority threshold,
    # same caveat as elsewhere in this function) pushed either edge
    # further than the unbaked reference's own detection sees.
    target_regions = []
    for region_start, region_end in _seed_regions:
        for r_start, r_end, _base, _length in find_homopolymer_runs("".join(reference_baked).upper(), 4):
            if r_start <= region_end and r_end >= region_start:
                region_start = min(region_start, r_start)
                region_end = max(region_end, r_end)
        target_regions.append((region_start, region_end))

    def in_target_regions(start, end):
        return any(start <= r_end and end >= r_start for r_start, r_end in target_regions)

    homopolymer_regions = [
        r for r in merge_runs(find_homopolymer_runs("".join(reference_baked).upper(), 4), 1)
        if in_target_regions(r[0], r[1])
    ]

    def region_containing(start, end):
        for r_start, r_end, _ in homopolymer_regions:
            if r_start <= start and end <= r_end:
                return (r_start, r_end)
        return None

    # marker_total_reads (computed above, excludes "Other sequences") is
    # also this reconciliation's own coverage denominator - a plain
    # Series, but its .get()/.items() work the same way a dict's would
    # for every use below.

    # Extract atomic deletion events directly from the original,
    # pre-explosion haplotype rows (see point 1 above).
    events = []
    for _, row in df_haplotypes.iterrows():
        tokens = str(row["sequence"]).split()
        del_positions = sorted(
            int(m.group(1)) for tok in tokens
            for m in [del_token_re.match(tok)] if m
        )
        run = []
        for p in del_positions:
            if run and p == run[-1] + 1:
                run.append(p)
            else:
                if run:
                    events.append({"positions": run, "marker": row["marker"], "total": row["total"]})
                run = [p]
        if run:
            events.append({"positions": run, "marker": row["marker"], "total": row["total"]})

    # Sum identical (marker, start, length) spans across different
    # original haplotype rows (e.g. the same 9bp deletion co-occurring
    # with several different, unrelated flanking SNPs).
    combined = {}
    for e in events:
        key = (e["marker"], e["positions"][0], len(e["positions"]))
        combined[key] = combined.get(key, 0) + e["total"]

    # Keep only events fully inside one reference homopolymer region -
    # that's where the repeat-shift math applies (see above).
    events = []
    for (marker, start, length), total in combined.items():
        region = region_containing(start, start + length - 1)
        if region is None:
            continue
        segment = "".join(reference_baked[start - 1: start - 1 + length])
        shifted_pos, shifted_segment = shift_deletion_right(reference_baked, start, segment)
        events.append({
            "marker": marker, "start": start, "length": length, "total": total,
            "region": region, "shifted_pos": shifted_pos, "shifted_segment": shifted_segment,
            "run_end": shifted_pos + len(shifted_segment) - 1,
        })

    lh_floor, lh_ceiling = lh_bounds_pct(length_heteroplasmy_threshold)

    merged_rows = []
    used_indices = set()

    # Boundary-run regions get their own independent, non-frame-baked
    # reporting (report_boundary_run above) for the RUN-LENGTH axis only.
    # Ordinary point substitutions strictly inside the extension run that
    # don't change length at all (e.g. T16189C, converting a T to a C
    # without adding or removing any bases) are NOT part of that axis -
    # they keep going through the general reconciliation below exactly
    # as before, including its own already-correct substitution-exclusion
    # handling (e.g. T16189C excluding T16189DEL's reads).
    events_for_combining = events
    # shared_frame: don't divert boundary-run regions away from the
    # general reconciliation at all - leave the region list empty so
    # every region below (16189-shaped or not) falls through to the same
    # single-frame handling 303-315-style plain homopolymer regions
    # already get.
    # Only chrM 16180 is boundary-run-shaped (a qualifying different-base
    # leading run in front of the extension run) - 303-315 has no such
    # leading run (see target_regions above), so it's never claimed here
    # and always falls through to the shared_frame-style handling below
    # regardless of mode, exactly like any other plain homopolymer run.
    boundary_run_regions = (
        [_target_br] if _target_br and homopolymer_reporting == "boundary_run" else []
    )
    # (marker, run_start, run_end) -> (dominant molecule's own label list,
    # pct) for every processed region, boundary-run-shaped or plain alike
    # - built up here and in the shared_frame/general path below, then
    # used in one final pass over `final` to append a "[dominant molecule
    # (...)]" tag to every row's own variant_note in that region (2026-
    # 09-22, user: "we think we can represent the information of the
    # dominant_only_shared frame in shared_frame option and dominant_
    # only_boundary_run in the boundary_run option ... let's include it
    # in the variant note, just the variant and percentage in brackets" -
    # then, after finding every row in a region carrying the SAME full
    # label list confusing (why doesn't -16193.1c also show -16193.1C
    # when that's exactly what it's reporting?): "I would like the
    # dominant molecule's variants to be matched to their rows: e.g. row
    # A16182a would only show A16182- (30.56%)". Keyed by (marker,
    # position) - extract_position(row's own sequence) against each
    # dominant row's own position, not the region's bounds, so a row only
    # ever gets the ONE dominant-molecule fact that's actually about the
    # same position as itself, matched regardless of case or major/minor
    # suffix (the row and the dominant call can genuinely differ there -
    # e.g. minor "A16182a" vs the dominant call's always-major "A16182-").
    # The full distribution stays exactly as it always was, row for row;
    # this is layered on top, not a replacement for it.
    dominant_summaries = {}

    for br in boundary_run_regions:
        leading_start, leading_end = br["leading"]
        extension_start, extension_end = br["extension"]
        extension_base = reference[extension_start - 1].upper()

        merged_rows.extend(report_boundary_run(
            final, df_haplotypes, marker_total_reads, reference,
            leading_start, leading_end, extension_start, extension_end,
            lh_floor, lh_ceiling
        ))
        dominant_rows = dominant_boundary_molecule(
            df_haplotypes, marker_total_reads, reference,
            leading_start, leading_end, extension_start, extension_end,
            lh_floor, lh_ceiling
        )
        # Keyed by the position extracted from the row's own sequence
        # label, not its "position" field - insertion rows across k=1..net
        # all share the SAME "position" (extension_end), only their own
        # decimal-anchored labels ("-16193.1C" vs "-16193.2C") actually
        # tell k=1 and k=2 apart.
        for r in dominant_rows:
            dominant_summaries[(r["marker"], extract_position(r["sequence"]))] = r

        # Drop only what report_boundary_run fully replaces: every raw
        # row at the leading run's own positions (sub or del alike - the
        # leading-run report above supersedes all of it), and the
        # extension run's own raw insertion tokens anchored at its far
        # end (the extension-run report's cumulative count supersedes
        # those specifically). Nothing else in 16184-16192 is touched.
        # Match plain integer-position substitution/deletion tokens
        # specifically - NOT via the generic extract_position helper,
        # which greedily parses "16182.1C"-style decimal insertion
        # tokens as if their integer part alone were a plain position.
        # Insertion tokens get their own, separate check just below:
        # anchored at the extension run's far end (superseded by its
        # cumulative count directly), OR anchored exactly at leading_end
        # (immediately adjacent to the extension run, folded into that
        # same cumulative count above for the common left=0 case - see
        # report_boundary_run's own left_edge comment for the general,
        # per-read-dynamic case this static check simplifies). An
        # insertion anchored any EARLIER in the leading run (e.g.
        # "16182.1C", with a reference A still sitting between it and
        # leading_end) is NOT part of either axis - a genuinely separate
        # event, left untouched here so it passes through and reports
        # itself at its own true position (2026-09-22, user: "the problem
        # here is that the insertion of 16182.1C occurs before the fourth
        # A and so we cannot shift it to 16193").
        plain_re = re.compile(r'^[ACGT](\d+)([ACGT]|DEL)$')
        for idx, row in final.iterrows():
            seq = str(row["sequence"])
            m = plain_re.match(seq)
            if m and leading_start <= int(m.group(1)) <= leading_end:
                used_indices.add(idx)
                continue
            ins_m = re.match(r'^(\d+)\.\d+([ACGT])$', seq)
            if ins_m:
                anchor, ins_base = int(ins_m.group(1)), ins_m.group(2)
                if anchor == extension_end:
                    used_indices.add(idx)
                elif anchor == leading_end and ins_base == extension_base:
                    used_indices.add(idx)

        # Suppress the general run-length combining loop below for
        # events report_boundary_run's extension-run count already fully
        # covers: one starting inside the leading run (e.g.
        # A16182DEL+A16183DEL - shift_deletion_right can't move an AA
        # deletion into the C-run, so it would otherwise get its own,
        # now-redundant combining-loop row right there in the leading
        # run), or one whose SHIFTED anchor lands inside the extension
        # run (e.g. T16189DEL, shifted to 16193 - already folded into
        # report_boundary_run's unified count). These events stay in the
        # full `events` list itself (used below) for two reasons: so
        # T16189C's own substitution-exclusion still correctly excludes
        # T16189DEL's reads, and so the del_mask step just below still
        # drops the now-superseded raw "T16189DEL"/"A16182DEL"/
        # "A16183DEL" single-position rows.
        events_for_combining = [
            e for e in events_for_combining
            if not (leading_start <= e["start"] <= leading_end)
            and not (extension_start <= e["run_end"] <= extension_end)
        ]

    events_by_run = {}
    for event in events_for_combining:
        events_by_run.setdefault((event["marker"], event["run_end"]), []).append(event)

    for (marker, run_end), run_events in events_by_run.items():
        coverage = marker_total_reads.get(marker, 0)
        max_len = max(len(e["shifted_segment"]) for e in run_events)
        # num_markers: how many distinct markers report this exact final
        # label elsewhere in `final` - default to 1 (this marker only)
        # since that's what these events are freshly derived from.
        for k in range(1, max_len + 1):
            contributing = [e for e in run_events if len(e["shifted_segment"]) >= k]
            reads = sum(e["total"] for e in contributing)
            if not reads or not coverage:
                continue
            pct = round(reads / coverage * 100, 1)
            if pct < lh_floor:
                continue
            pos = run_end - k + 1
            base = reference_baked[pos - 1]
            is_major = pct >= lh_ceiling
            # Explanatory variant_note suppressed - not deleted, same as
            # report_boundary_run's own notes above (2026-09-22, user:
            # "too much information ... I don't want to see it").
            # detail = ", ".join(
            #     f"{e['total']} read{'s' if e['total'] != 1 else ''} with a deletion at position"
            #     + (f"s {e['start']}-{e['start'] + e['length'] - 1}" if e["length"] > 1 else f" {e['start']}")
            #     for e in contributing
            # )
            # note = f"{reads} of {coverage} reads ({pct}%) are at least {k} base{'s' if k != 1 else ''} shorter than the reference here"
            # note += f", combining {len(contributing)} separate deletions ({detail})" if len(contributing) > 1 else f" ({detail})"
            merged_rows.append({
                "sequence": f"{base}{pos}" + ("-" if is_major else base.lower()),
                "total": reads,
                "interpolated_total_coverage": coverage,
                "is_noise_or_low_frq": False,
                "num_markers": 1,
                "variant_frequency": pct,
                "marker": marker,
                "position": float(pos),
            })

    # --- Insertion-side mirror of the deletion combining loop above ---
    #
    # Deletions get real reconciliation (atomic event extraction, repeat-
    # shift, cumulative combining); insertions in this general path never
    # did - a raw FDSTOOLS insertion token (e.g. "-16193.1c") just passed
    # through unreconciled and unexplained, silently, since nothing here
    # ever looked at it. That's the actual gap behind homopolymer_
    # reporting=shared_frame's missing variant_note, not something the
    # boundary-run-specific fix above touches (shared_frame skips that
    # entirely). This closes it the same way deletions are already
    # closed: per individual homopolymer run (not the merged, possibly
    # mixed-base regions above - an insertion's own base has to match the
    # ONE run it's extending, so this needs single-base runs, decomposing
    # a boundary-run-shaped region into its separate leading/extension
    # pieces the same way report_boundary_run's own two axes do), each
    # read's net extra-base contribution is: insertions of the run's own
    # base anchored anywhere inside it, plus a walk past the run's own
    # far end for substitutions converting to that base (mirrors report_
    # boundary_run's right_sub - the same real case, HG02389's "A16194C
    # 16194.1C", applies here too whenever shared_frame is active), plus
    # any insertion anchored at wherever that walk lands.
    #
    # Boundary-run-shaped regions are skipped here when they were already
    # claimed by report_boundary_run above (i.e. whenever homopolymer_
    # reporting is boundary_run) - nothing to reconcile twice. For
    # shared_frame, boundary_run_regions is empty, so nothing is excluded,
    # and these get exactly the same treatment as any other homopolymer
    # run.
    ins_re_general = re.compile(r'^(\d+)\.(\d+)([ACGT])$')
    claimed_spans = [span for br in boundary_run_regions for span in (br["leading"], br["extension"])]

    def is_claimed(start, end):
        # Overlap, not containment: baking can extend a run's own boundary
        # (e.g. a homoplasmic A16183C bakes position 16183 to C, merging
        # what find_homopolymer_runs sees into one 16183-16193 run) so
        # that it straddles both the leading and extension spans without
        # being fully contained in either - still needs excluding, since
        # report_boundary_run already covers that whole stretch correctly
        # on its own, per-read terms.
        return any(start <= c_end and end >= c_start for c_start, c_end in claimed_spans)

    # Restricted to target_regions (16180 and 303-315) same as
    # boundary_run_regions/homopolymer_regions above - every other plain
    # 4+ homopolymer run genome-wide is left untouched, passed through
    # exactly as FDSTOOLS reported it.
    homopolymer_runs = [
        (r_start, r_end, base, length)
        for r_start, r_end, base, length in find_homopolymer_runs("".join(reference_baked).upper(), 4)
        if not is_claimed(r_start, r_end) and in_target_regions(r_start, r_end)
    ]

    ins_events_by_run = {}
    ins_used_indices = set()
    for _, row in df_haplotypes.iterrows():
        seq_str = str(row["sequence"]).strip()
        if seq_str == "Other sequences":
            continue
        tokens = seq_str.split()
        marker, n = row["marker"], row["total"]
        subbed_row = {int(m.group(2)): m.group(3) for tok in tokens for m in [sub_pattern.match(tok)] if m}

        for run_start, run_end, run_base, _length in homopolymer_runs:
            ins_by_anchor = {}
            for tok in tokens:
                m = ins_re_general.match(tok)
                if m and m.group(3) == run_base and run_start <= int(m.group(1)) <= run_end:
                    a = int(m.group(1))
                    ins_by_anchor[a] = max(ins_by_anchor.get(a, 0), int(m.group(2)))
            inside = sum(ins_by_anchor.values())

            # Mirror pair, extending the run outward on either side for
            # this one read: a substitution immediately adjacent to the
            # run's own start/end, converting to the run's base, plus any
            # insertion anchored at wherever that walk lands. Symmetric
            # with report_boundary_run's leading-run "left" term and this
            # module's own right_sub above - a per-read, non-baked
            # substitution at either edge (e.g. a minority A16182C next
            # to a baked-in A16183C) is otherwise invisible here, since
            # baking only extends the run's REFERENCE-DERIVED boundary
            # for substitutions common enough to be baked at all.
            left_sub = 0
            p = run_start - 1
            while subbed_row.get(p) == run_base:
                left_sub += 1
                p -= 1
            left_edge = run_start - left_sub

            right_sub = 0
            p = run_end + 1
            while subbed_row.get(p) == run_base:
                right_sub += 1
                p += 1
            right_edge = run_end + right_sub

            # Only counted when the walk actually left the run (left_sub/
            # right_sub > 0): when it doesn't, left_edge-1 or right_edge
            # collapses onto run_start-1 (always outside `inside`'s own
            # range, harmless) or exactly run_end (NOT outside it - a
            # plain "16193.1C" with no edge substitution would otherwise
            # be counted once by `inside` and a second time here). Each
            # side maxed independently first (FDSTOOLS' own multi-token
            # convention for "K inserted bases", e.g. "16193.1C 16193.2C
            # 16193.3C" meaning 3, not 1+2+3), then the two sides summed,
            # since a real left-side and right-side extension on the same
            # read are genuinely independent events.
            left_edge_ins = 0
            if left_sub > 0:
                for tok in tokens:
                    m = ins_re_general.match(tok)
                    if m and int(m.group(1)) == left_edge - 1:
                        left_edge_ins = max(left_edge_ins, int(m.group(2)))
            right_edge_ins = 0
            if right_sub > 0:
                for tok in tokens:
                    m = ins_re_general.match(tok)
                    if m and int(m.group(1)) == right_edge:
                        right_edge_ins = max(right_edge_ins, int(m.group(2)))

            net_ins = inside + left_sub + right_sub + left_edge_ins + right_edge_ins
            if net_ins > 0:
                ins_events_by_run.setdefault((marker, run_start, run_end, run_base), []).append((net_ins, n))

    for (marker, run_start, run_end, run_base), entries in ins_events_by_run.items():
        coverage = marker_total_reads.get(marker, 0)
        if not coverage:
            continue
        max_net = max(net for net, _n in entries)
        for k in range(1, max_net + 1):
            reads = sum(n for net, n in entries if net >= k)
            if not reads:
                continue
            pct = round(reads / coverage * 100, 1)
            if pct < lh_floor:
                continue
            is_major = pct >= lh_ceiling
            base_label = run_base if is_major else run_base.lower()
            # Explanatory variant_note suppressed - not deleted, same as
            # the other combining loops above.
            # note = (f"{reads} of {coverage} reads ({pct}%) have at least {k} extra "
            #         f"{run_base}{'s' if k != 1 else ''} in the {run_start}-{run_end} run - "
            #         f"from insertions anchored inside the run, substitutions extending past "
            #         f"its own far end, or both combined")
            merged_rows.append({
                "sequence": f"-{run_end}.{k}{base_label}",
                "total": reads,
                "interpolated_total_coverage": coverage,
                "is_noise_or_low_frq": False,
                "num_markers": 1,
                "variant_frequency": pct,
                "marker": marker,
                "position": float(run_end),
            })

    for run_start, run_end, run_base, _length in homopolymer_runs:
        # Drop the raw insertion rows this now supersedes: any exploded
        # "POS.KBASE" row anchored inside this run (or at the exact right
        # edge some read's own right_sub walk could reach - matched
        # loosely here since `final` no longer carries per-read context,
        # just anchored at run_end itself, the common case). Matched on
        # position and base alone, not marker - a raw row only exists for
        # a marker whose own reads actually produced it, so it's already
        # implicitly that marker's own row regardless.
        for idx, frow in final.iterrows():
            fm = re.match(r'^(\d+)\.\d+([ACGT])$', str(frow["sequence"]))
            if fm and fm.group(2) == run_base and run_start <= int(fm.group(1)) <= run_end:
                ins_used_indices.add(idx)

    used_indices |= ins_used_indices

    # Drop the exploded DEL rows now superseded by the combined events
    # above - i.e. any DEL row, for a marker with at least one processed
    # event, whose position falls in that same homopolymer region.
    processed_regions = {(e["marker"], e["region"]) for e in events}
    del_mask = final["sequence"].astype(str).str.match(r'^[ACGT]\d+DEL$')
    if processed_regions:
        for idx, row in final[del_mask].iterrows():
            pos = int(re.search(r'\d+', row["sequence"]).group())
            for marker, (r_start, r_end) in processed_regions:
                if row["marker"] == marker and r_start <= pos <= r_end:
                    used_indices.add(idx)
                    break

    # shared_frame's own dominant-molecule summary, computed jointly
    # (insertion and deletion together per read) exactly as before -
    # folded into dominant_summaries for the tagging pass below instead
    # of replacing any of the full-distribution rows above. Harmless
    # no-op under boundary_run mode: homopolymer_runs is empty there
    # (this region is_claimed by boundary_run_regions instead).
    baked_positions_by_run = find_baked_positions_by_run(reference, homopolymer_runs)
    dominant_rows = dominant_shared_frame_molecule(
        df_haplotypes, marker_total_reads, reference, reference_baked,
        homopolymer_runs, baked_positions_by_run, lh_floor, lh_ceiling
    )
    # Keyed by the position extracted from the row's own sequence label,
    # same reasoning as the boundary_run-shaped dispatch above.
    for r in dominant_rows:
        dominant_summaries[(r["marker"], extract_position(r["sequence"]))] = r

    # Substitution rows whose position falls in a processed homopolymer
    # region: the denominator always stays the marker's full coverage
    # (not shrunk by the co-located deletion's reads), so every row in
    # the amplicon shares one consistent total. Whether the deletion's
    # reads also join the NUMERATOR depends on whether any reads at all
    # still retain the true reference base here: if not - every read
    # either substitutes or deletes, e.g. T16189C's 77 + T16189DEL's 29
    # = coverage exactly - the deletion reads are just as much "not
    # reference" evidence as the substitution reads, so folding them in
    # correctly reports "100%, nothing here matches reference" instead
    # of artificially withholding them. If some reads DO still show true
    # reference (e.g. A16182's 64 + A16182DEL's 8 = 72, short of
    # coverage 106 by the 34 reads that keep reference A), the deletion
    # reads are a genuine third category and stay out of the numerator,
    # so the real ref/alt split isn't distorted - they still count in
    # the denominator, though, since that's now always the full coverage.
    if "variant_note" not in final.columns:
        final["variant_note"] = pd.NA
    sub_mask = final["sequence"].astype(str).str.match(r'^[ACGT]\d+[ACGT]$')
    for idx, row in final[sub_mask].iterrows():
        m = sub_pattern.match(row["sequence"])
        pos = int(m.group(2))
        covering = [e for e in events if e["marker"] == row["marker"] and e["start"] <= pos <= e["start"] + e["length"] - 1]
        if not covering:
            continue
        del_total = sum(e["total"] for e in covering)
        coverage = marker_total_reads.get(row["marker"])
        if not coverage:
            continue
        alt_total = final.at[idx, "total"]
        ref_total = coverage - alt_total - del_total
        fold_in = ref_total <= 0
        numerator = alt_total + del_total if fold_in else alt_total
        new_freq = round(numerator / coverage * 100, 2)
        # Numbers only, no explanatory variant_note text - the FMP label,
        # total, and frequency still fully reflect the fold-in logic
        # above, just without spelling out why in prose. Suppressed, not
        # deleted (2026-09-22, user: "too much to follow" / "let's not
        # delete it" / "suppress it") - kept below so it can be turned
        # back on.
        # del_desc = ", ".join(
        #     f"{e['total']} read{'s' if e['total'] != 1 else ''} with "
        #     + "+".join(f"{reference[p - 1]}{p}DEL" for p in range(e["start"], e["start"] + e["length"]))
        #     for e in covering
        # )
        # note = (
        #     f"{del_total} of {coverage} reads ({round(del_total / coverage * 100, 1)}%) have a deletion "
        #     f"at this position instead ({del_desc}) - counted here too, since no reads at this position "
        #     f"still show the reference base"
        #     if fold_in else
        #     f"{del_total} of {coverage} reads ({round(del_total / coverage * 100, 1)}%) have a deletion "
        #     f"at this position instead ({del_desc}) - not counted toward this percentage, "
        #     f"but still part of the {coverage} total"
        # )
        final.at[idx, "total"] = numerator
        final.at[idx, "interpolated_total_coverage"] = coverage
        final.at[idx, "variant_frequency"] = new_freq

    # --- Fallback: original single-pair naive-position merge, for
    # whatever del+sub pairs weren't handled by the homopolymer-region
    # logic above (e.g. the 8281-8289 tandem-repeat deletion). ---
    remaining = final.drop(index=used_indices)
    for pos, group in remaining.groupby("position"):
        if group.shape[0] != 2:
            continue
        del_row = group[group["sequence"].str.endswith("DEL")]
        sub_row = group[~group["sequence"].str.endswith("DEL")]
        if del_row.empty or sub_row.empty:
            continue
        del_idx, sub_idx = del_row.index[0], sub_row.index[0]
        if del_idx in used_indices or sub_idx in used_indices:
            continue
        match = re.match(r'([ACGT])(\d+)([ACGT])', sub_row["sequence"].iloc[0])
        if not match:
            print(f"WARNING: Could not parse sequence: {sub_row['sequence'].iloc[0]}")
            continue
        ref, pos_str, alt = match.groups()
        sub_pos = int(pos_str)

        sub_total = sub_row["total"].iloc[0]
        del_total = del_row["total"].iloc[0]
        coverage = del_row["interpolated_total_coverage"].iloc[0]

        del_freq = round(del_total / coverage * 100, 1) if coverage else 0
        del_seq_name = del_row['sequence'].iloc[0]
        deletion_reportable = del_freq >= lh_floor

        if deletion_reportable:
            reported_total = sub_total + del_total
            sub_note = (f"includes {del_total} of {coverage} reads ({del_freq}%) with {del_seq_name} "
                        f"- also reported as its own row below")
        else:
            reported_total = sub_total
            sub_note = (f"excludes {del_total} of {coverage} reads ({del_freq}%) with {del_seq_name} "
                        f"- below the {lh_floor:g}% reporting threshold, so not counted here or shown separately")
        reported_freq = round(reported_total / coverage * 100, 1) if coverage else 0

        num_markers = sub_row['num_markers'].iloc[0]
        marker = sub_row["marker"].iloc[0]

        merged_rows.append({
            "sequence": f"{ref}{pos_str}{alt}",
            "total": reported_total,
            "interpolated_total_coverage": coverage,
            "is_noise_or_low_frq": False,
            "num_markers": num_markers,
            "variant_frequency": reported_freq,
            "variant_note": sub_note,
            "marker": marker,
            "position": float(pos)
        })

        reference_with_sub = list(reference_baked)
        del_pos, del_segment = shift_deletion_right(reference_with_sub, sub_pos, reference_with_sub[sub_pos - 1])

        if deletion_reportable:
            is_major = del_freq >= lh_ceiling
            del_seq = f"{del_segment}{del_pos}" + ("-" if is_major else del_segment.lower())
            merged_rows.append({
                "sequence": del_seq,
                "total": del_total,
                "interpolated_total_coverage": coverage,
                "is_noise_or_low_frq": False,
                "num_markers": num_markers,
                "variant_frequency": del_freq,
                "variant_note": (f"{del_total}/{coverage} ({del_freq}%) of {ref}{pos_str}{alt} reads - "
                                 f"deletion shifted from naive position {pos_str}"),
                "marker": marker,
                "position": float(del_pos)
            })

        used_indices.update([del_idx, sub_idx])

    # Drop merged original rows from final DataFrame
    final = final.drop(index=used_indices)
    if merged_rows:
        final = pd.concat([final, pd.DataFrame(merged_rows)], ignore_index=True)

    # Tag each row with just its OWN matching dominant-molecule fact, not
    # the whole region's label list - matched by the position extracted
    # from each side's own sequence label (case- and major/minor-suffix-
    # insensitive, e.g. minor "A16182a" still matches the dominant call's
    # always-major "A16182-"; case-insensitive so "-16193.1c" matches the
    # dominant call's own "-16193.1C"), not the region's bounds (2026-09-
    # 22, user: "I'm confused why we wouldn't also include in the rows
    # 16193.1c and 16193.2c since the dominant molecule contains -16193.1C
    # and -16193.2C ... instead of showing the same message I would like
    # the dominant molecule's variants to be matched to their rows: e.g.
    # row A16182a would only show A16182- (30.56%)"). A row with no
    # matching dominant-molecule position (e.g. -16193.3c when the
    # dominant molecule's own net only reaches 2) gets nothing.
    if "variant_note" not in final.columns:
        final["variant_note"] = pd.NA
    if dominant_summaries:
        matched_keys = set()
        for idx, row in final.iterrows():
            key = (row["marker"], extract_position(str(row["sequence"])))
            match = dominant_summaries.get(key)
            if not match:
                continue
            matched_keys.add(key)
            tag = f"{match['sequence']} ({match['variant_frequency']}%)"
            existing = final.at[idx, "variant_note"]
            if pd.isna(existing) or not str(existing).strip():
                final.at[idx, "variant_note"] = tag
            else:
                final.at[idx, "variant_note"] = f"{existing} {tag}"

        # A dominant-molecule fact with no existing row to attach to
        # (e.g. HG01799's "C16193-": report_boundary_run's own full-
        # distribution model has no notation at all for a run coming out
        # net SHORTER than reference, only "at least k extra") would
        # otherwise be silently computed and then lost - add it as its
        # own row instead of dropping a real finding just because the
        # full-distribution mechanism happens to have a gap for that
        # particular shape.
        missing = [r for key, r in dominant_summaries.items() if key not in matched_keys]
        if missing:
            # A new row becomes a normal, standalone FMP entry sitting
            # next to every other row in the file - it needs their same
            # major/minor convention (uppercase "-" at/above lh_ceiling,
            # lowercase reference base below it), not the dominant
            # functions' own deliberately-always-major convention (which
            # is right for a tag describing "this IS the dominant
            # molecule's genotype", but wrong for a label a reader
            # compares against every other row's own frequency). General
            # over any position/base/sample - driven only by lh_ceiling
            # and the row's own frequency, not anything specific to this
            # one case (2026-09-22, user: "are you hardcoding this? we
            # want a solution that works generally"). Every "-"-shaped
            # row any of these dominant functions can produce (leading-
            # run loss, either one's own net<0 branch) is this same
            # plain "{base}{pos}-" shape, so one general regex covers all
            # of them.
            plain_del_re = re.compile(r'^([ACGT])(\d+)-$')
            new_rows = []
            for r in missing:
                seq = r["sequence"]
                m = plain_del_re.match(seq)
                if m and r["variant_frequency"] < lh_ceiling:
                    base, pos = m.group(1), m.group(2)
                    seq = f"{base}{pos}{base.lower()}"
                new_rows.append({**r, "sequence": seq, "variant_note": f"{r['sequence']} ({r['variant_frequency']}%)"})
            final = pd.concat([final, pd.DataFrame(new_rows)], ignore_index=True)

    # If nothing remains after filtering (e.g., H2O / No data), write empty output and stop
    if final.empty:
        pd.DataFrame(columns=[
            "FDSTOOLS", "vf_FDS", "rd_FDS", "interpolated_total_coverage",
            "variant_note", "marker", "marker_range", "num_markers"
        ]).to_csv(output_file, sep="\t", index=False)
        print(f"Output written to {output_file} (no variants)")
        return
    resolved = final.apply(
        lambda row: resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold, IUPAC_CODES),
        axis=1
    )

    # Force to a plain 1D Series of strings
    final["sequence"] = pd.Series(resolved, index=final.index).astype(str)

    # Drop length variants below the symmetric length-heteroplasmy floor entirely
    final = final[final["sequence"] != LH_DROP_SENTINEL]
    
    # final["sequence"] = final.apply(
    #     lambda row: resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold, IUPAC_CODES),
    #     axis=1
    # )

    final = pd.concat([final, single_low_coverage], ignore_index=False)

    final["marker_range"] = final["marker"].map(marker_map)
    final["position"] = final["marker"].apply(extract_position)
    final["position"] = final["position"].apply(adjust_circular_position)
    final = final.sort_values(by="position").drop(columns=["position"])

    # Rename columns
    final = final.rename(columns={
        "sequence": "FDSTOOLS",
        "variant_frequency": "vf_FDS",
        "total": "rd_FDS"
    })

    # Reorder columns
    desired_order = [
        "FDSTOOLS",
        "vf_FDS",
        "rd_FDS",
        "interpolated_total_coverage",
        "variant_note",
        "marker",
        "marker_range",
        "num_markers"
    ]
    existing_columns = [col for col in desired_order if col in final.columns]
    final = final[existing_columns + [col for col in final.columns if col not in existing_columns]]

    # if final is None or final.empty:
    #     pd.DataFrame(columns=[
    #         "FDSTOOLS", "vf_FDS", "rd_FDS", "interpolated_total_coverage",
    #         "variant_note", "marker", "marker_range", "num_markers"
    #     ]).to_csv(output_file, sep="\t", index=False)
    #     print(f"Output written to {output_file} (no variants)")
    #     return
    final.to_csv(output_file, sep="\t", index=False)
    print(f"Output written to {output_file}")



# CLI wrapper
def main():
    parser = argparse.ArgumentParser(description="Process FDSTools SAST TSV file with heteroplasmy handling.")
    parser.add_argument("input", help="Input SAST TSV file")
    parser.add_argument("output", help="Output processed file")
    parser.add_argument("--marker_map", help="Path to marker map file")
    parser.add_argument("--reference", required=True, help="Reference genome in FASTA format, used to repeat-shift deletions separated out of a substitution+deletion pairing")
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minimum variant frequency threshold")
    parser.add_argument("--depth", type=int, default=10, help="Read depth threshold for low coverage")
    parser.add_argument("--lh_thresh", type=float, default=10.0, help="Symmetric length-heteroplasmy threshold (floor and, via 100-threshold, ceiling), e.g. 10.0 -> report only 10-90%%, lowercase in between, major above 90%%")
    parser.add_argument("--homopolymer_reporting",
                         choices=["shared_frame", "boundary_run"],
                         default="boundary_run",
                         help="Which underlying model reports homopolymer run lengths, restricted to just chrM 16180 through "
                              "the end of its poly-C tract (16189 T>C's own consequence) - the only region any of this has "
                              "been validated against. Every other homopolymer run genome-wide (310, the other 8 boundary-run "
                              "regions, any other plain 4+ run) is passed through exactly as FDSTOOLS reported it. "
                              "shared_frame bakes everything into one borrowed reference frame. boundary_run (default) splits "
                              "the region's leading run and different-base extension run, 16180-16183 / 16184-16193, into two "
                              "independent axes instead. Either way, every row in the region also gets the single most common "
                              "molecule's own genotype appended to its variant_note as a bracketed tag (e.g. \"[dominant "
                              "molecule (30.56%%): A16182C, A16183C, T16189C, -16193.1C, -16193.2C]\"), per W. Parson's "
                              "reading of the ISFG convention - the full distribution stays exactly as it always was, this is "
                              "layered on top, not a replacement for it.")
    args = parser.parse_args()

    try:
        process_fdstools_sast(
            file_path=args.input,
            marker_map_path=args.marker_map,
            output_file=args.output,
            reference_fasta=args.reference,
            min_variant_frequency_pct=args.min_vf,
            depth_threshold=args.depth,
            length_heteroplasmy_threshold=args.lh_thresh,
            homopolymer_reporting=args.homopolymer_reporting
        )
    except Exception:
        traceback.print_exc()
        sys.exit(1)
    # except Exception as e:
    #     print(f"Error: {e}", file=sys.stderr)
    #     sys.exit(1)

if __name__ == "__main__":
    main()
