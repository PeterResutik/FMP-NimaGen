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

    if 'DEL' in seq:
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
            if loss_pct >= lh_ceiling:
                rows.append({
                    "sequence": f"{ref_here}{pos}-",
                    "total": coverage - ref_total,
                    "interpolated_total_coverage": coverage,
                    "is_noise_or_low_frq": False,
                    "num_markers": 1,
                    "variant_frequency": loss_pct,
                    "variant_note": f"reference {ref_here} observed at only {int(ref_total)}/{coverage} reads here "
                                     f"- the rest are reported as part of the {extension_start}-{extension_end} run below",
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
                    "variant_note": f"reference {ref_here} retained at {int(ref_total)}/{coverage} reads here "
                                     f"- the rest are reported as part of the {extension_start}-{extension_end} run below",
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

        deleted = {int(m.group(1)) for tok in tokens for m in [del_re.match(tok)] if m}
        internal = sum(1 for d in deleted if extension_start <= d <= extension_end)

        # An insertion of the extension run's own base anchored ANYWHERE
        # in the leading run (e.g. "16182.1C") extends the same run just
        # as much as one anchored at its far end (e.g. "16193.1C") - both
        # are "one more C than reference has", full stop, regardless of
        # which position FDSTOOLS' own per-haplotype notation happened to
        # anchor it at. Folding it in here is also what makes a
        # cancelling insertion+deletion pair (HG01799's "16182.1C" +
        # "T16189DEL" together) correctly net to 0 - a normal-length run
        # - instead of registering as a pure -1 loss.
        leading_insertion_by_anchor = {}
        for tok in tokens:
            m = ins_re.match(tok)
            if m and m.group(3) == extension_base and leading_start <= int(m.group(1)) <= leading_end:
                anchor = int(m.group(1))
                leading_insertion_by_anchor[anchor] = max(leading_insertion_by_anchor.get(anchor, 0), int(m.group(2)))
        leading_insertion = sum(leading_insertion_by_anchor.values())

        right = 0
        for tok in tokens:
            m = ins_re.match(tok)
            if m and int(m.group(1)) == extension_end:
                right = max(right, int(m.group(2)))

        net = left + leading_insertion - internal + right
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
            rows.append({
                "sequence": f"-{extension_end}.{k}{base}",
                "total": reads,
                "interpolated_total_coverage": coverage,
                "is_noise_or_low_frq": False,
                "num_markers": 1,
                "variant_frequency": pct,
                "variant_note": (f"{reads} of {coverage} reads ({pct}%) have at least {k} extra "
                                  f"{extension_base}{'s' if k != 1 else ''} in the {extension_start}-{extension_end} "
                                  f"run - from substitutions at the run's own leading edge "
                                  f"({leading_start}-{leading_end}), insertions past {extension_end}, or both combined"),
                "marker": marker,
                "position": float(extension_end),
            })
            k += 1

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
def process_fdstools_sast(file_path, marker_map_path, output_file, reference_fasta, min_variant_frequency_pct=5.0, depth_threshold=10, length_heteroplasmy_threshold=90.0):
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
    # df["is_noise_or_low_frq"] = df["sequence"].isin(["Other sequences"]) | (df["total_mp_sum"] < min_variant_frequency_pct)
    # clean_total_per_marker = df[~df["is_noise_or_low_frq"]].groupby("marker")["total"].sum().rename("total_wo_noise_or_low_frq")
    # df = df.merge(clean_total_per_marker, on="marker", how="left")
    # df["variant_frequency_wo_noise_or_low_frq"] = (df["total"] / df["total_wo_noise_or_low_frq"] * 100).round(2)

    # df = df.assign(sequence=df["sequence"].str.split()).explode("sequence").reset_index(drop=True)
    # drop_seqs = ["Other", "sequences", "REF", "N3107DEL"]
    # df = df[(~df["sequence"].isin(drop_seqs)) & (df["total_mp_sum"] >= min_variant_frequency_pct)].copy()

    # Step 5: Split multiple variants
    df = df.assign(sequence=df["sequence"].str.split())
    df = df.explode("sequence").reset_index(drop=True)
    
    # df["interpolated_total_coverage"] = (np.ceil(df["total"] / (df["total_mp_sum"] / 100))).astype("Int64")

    denom = (df["total_mp_sum"] / 100).replace(0, np.nan)  # avoid division by zero
    interp = np.ceil(df["total"] / denom)  # will be NaN where denom was 0
    df["interpolated_total_coverage"] = pd.to_numeric(interp, errors="coerce").fillna(0).astype("Int64")

    grouped = df.groupby(["marker", "sequence"], as_index=False).agg(
        total=("total", "sum"),
        total_mp_sum=("total_mp_sum", "sum"),
        interpolated_total_coverage=("interpolated_total_coverage", "max"),
        # is_noise_or_low_frq=("is_noise_or_low_frq", "first"),
        # total_wo_noise_or_low_frq=("total_wo_noise_or_low_frq", "first"),
        # variant_frequency_wo_noise_or_low_frq=("variant_frequency_wo_noise_or_low_frq", "sum")
    )

    # Extract "Other" sequence coverage per marker
    other_per_marker = grouped[grouped["sequence"] == "Other"][["marker", "total"]]
    other_per_marker = other_per_marker.rename(columns={"total": "other_coverage"})

    # Merge with grouped data
    grouped = grouped.merge(other_per_marker, on="marker", how="left")
    grouped["other_coverage"] = grouped["other_coverage"].fillna(0)

    # Subtract "Other" sequence coverage from interpolated_total_coverage
    grouped["adjusted_coverage"] = grouped["interpolated_total_coverage"] - grouped["other_coverage"]

    # Ensure adjusted coverage is not negative or zero (to avoid division by zero)
    grouped["interpolated_total_coverage"] = grouped["adjusted_coverage"].clip(lower=1)


    final = grouped.groupby("sequence", as_index=False).agg(
        marker=("marker", "first"),
        total=("total", "sum"),
        interpolated_total_coverage=("interpolated_total_coverage", "sum"),
        # is_noise_or_low_frq=("is_noise_or_low_frq", "first"),
        # total_wo_noise_or_low_frq=("total_wo_noise_or_low_frq", "sum"),
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

    homopolymer_regions = merge_runs(find_homopolymer_runs("".join(reference_baked).upper(), 4), 1)

    def region_containing(start, end):
        for r_start, r_end, _ in homopolymer_regions:
            if r_start <= start and end <= r_end:
                return (r_start, r_end)
        return None

    # True per-marker read total (excludes the "Other sequences" catch-
    # all bucket) - used as the reconciliation's coverage denominator
    # instead of each position's own interpolated coverage, which can
    # wobble by rounding error on FDSTOOLS' own 3-sig-fig percentages.
    is_other = df_haplotypes["sequence"].astype(str).str.strip() == "Other sequences"
    marker_total_reads = df_haplotypes[~is_other].groupby("marker")["total"].sum().to_dict()

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
    boundary_run_regions = find_boundary_run_regions("".join(reference).upper())
    for br in boundary_run_regions:
        leading_start, leading_end = br["leading"]
        extension_start, extension_end = br["extension"]
        extension_base = reference[extension_start - 1].upper()

        merged_rows.extend(report_boundary_run(
            final, df_haplotypes, marker_total_reads, reference,
            leading_start, leading_end, extension_start, extension_end,
            lh_floor, lh_ceiling
        ))

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
        # cumulative count directly), OR anchored in the leading run but
        # inserting the EXTENSION run's own base (e.g. "16182.1C" -
        # folded into that same cumulative count above, so it must not
        # ALSO still appear here as its own separate row). An insertion
        # of some OTHER base in the leading run isn't part of either
        # axis this function reports and is left untouched.
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
                elif leading_start <= anchor <= leading_end and ins_base == extension_base:
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
            detail = ", ".join(
                f"{e['total']} read{'s' if e['total'] != 1 else ''} with a deletion at position"
                + (f"s {e['start']}-{e['start'] + e['length'] - 1}" if e["length"] > 1 else f" {e['start']}")
                for e in contributing
            )
            note = f"{reads} of {coverage} reads ({pct}%) are at least {k} base{'s' if k != 1 else ''} shorter than the reference here"
            note += f", combining {len(contributing)} separate deletions ({detail})" if len(contributing) > 1 else f" ({detail})"
            merged_rows.append({
                "sequence": f"{base}{pos}" + ("-" if is_major else base.lower()),
                "total": reads,
                "interpolated_total_coverage": coverage,
                "is_noise_or_low_frq": False,
                "num_markers": 1,
                "variant_frequency": pct,
                "variant_note": note,
                "marker": marker,
                "position": float(pos),
            })

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
        del_desc = ", ".join(
            f"{e['total']} read{'s' if e['total'] != 1 else ''} with "
            + "+".join(f"{reference[p - 1]}{p}DEL" for p in range(e["start"], e["start"] + e["length"]))
            for e in covering
        )
        final.at[idx, "total"] = numerator
        final.at[idx, "interpolated_total_coverage"] = coverage
        final.at[idx, "variant_frequency"] = new_freq
        if fold_in:
            final.at[idx, "variant_note"] = (
                f"{del_total} of {coverage} reads ({round(del_total / coverage * 100, 1)}%) have a deletion "
                f"at this position instead ({del_desc}) - counted here too, since no reads at this position "
                f"still show the reference base"
            )
        else:
            final.at[idx, "variant_note"] = (
                f"{del_total} of {coverage} reads ({round(del_total / coverage * 100, 1)}%) have a deletion "
                f"at this position instead ({del_desc}) - not counted toward this percentage, "
                f"but still part of the {coverage} total"
            )

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
    args = parser.parse_args()

    try:
        process_fdstools_sast(
            file_path=args.input,
            marker_map_path=args.marker_map,
            output_file=args.output,
            reference_fasta=args.reference,
            min_variant_frequency_pct=args.min_vf,
            depth_threshold=args.depth,
            length_heteroplasmy_threshold=args.lh_thresh
        )
    except Exception:
        traceback.print_exc()
        sys.exit(1)
    # except Exception as e:
    #     print(f"Error: {e}", file=sys.stderr)
    #     sys.exit(1)

if __name__ == "__main__":
    main()
