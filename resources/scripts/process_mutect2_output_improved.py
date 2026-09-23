#!/usr/bin/env python3

import pandas as pd
import re
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from call_repeat_regions import MUTECT2_TARGET_REGIONS
from reference_utils import (
    load_reference_list, rightmost_repeat_position, shift_insertion_right,
    shift_deletion_right,
)
from heteroplasmy_thresholds import lh_bounds

IUPAC_CODES = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["A", "C"]): "M",
    frozenset(["G", "T"]): "K",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W"
}

# rCRS_NimaGen.fasta is the real, 16569bp circular rCRS with a 54bp copy of
# chrM:1-54 appended, so the last (origin-spanning) amplicon can map
# linearly instead of soft-clipping across the wrap. That appended stretch
# isn't a real extra locus - a Mutect2 POS in 16570-16623 IS the same
# physical base as POS-16569, just observed through the wraparound
# amplicon's copy of the reference. Every position that reaches a
# human-readable label must be converted back before it does, or it comes
# out as an invalid >16569 coordinate that can't match FDSTOOLS' own
# (correctly wrapped) labels in the downstream merge step.
TRUE_MT_LENGTH = 16569

def wrap_circular_position(pos, true_length=TRUE_MT_LENGTH):
    return pos - true_length if pos > true_length else pos

def apply_snp(pos, ref, var, var_level, reference, min_variant_frequency):
    is_major = var_level >= 1 - min_variant_frequency
    formatted = []
    for i, (r, v) in enumerate(zip(ref, var)):
        sub_pos = wrap_circular_position(pos + i)
        if is_major:
            formatted.append(f"{r}{sub_pos}{v}")
        else:
            code = IUPAC_CODES.get(frozenset([r, v]), f"{r}/{v}")
            formatted.append(f"{r}{sub_pos}{code}")

    return " ".join(formatted), ("SNP" if is_major else "PHP")

def apply_insertion(pos, ref, var, var_level, reference, length_heteroplasmy_threshold):
    floor, ceiling = lh_bounds(length_heteroplasmy_threshold)
    if var_level < floor:
        return None, "BELOW_LH_FLOOR"
    inserted_segment = var[len(ref):]
    pos, segment = shift_insertion_right(reference, pos, inserted_segment)
    pos = wrap_circular_position(pos)
    is_major = var_level >= ceiling
    variant_parts = [
        f"-{pos}.{i+1}{(b if is_major else b.lower())}"
        for i, b in enumerate(segment)
    ]
    updated_type = "INS" if is_major else "LHP"
    return " ".join(variant_parts), updated_type

def apply_deletion(pos, ref, var, var_level, reference, length_heteroplasmy_threshold):
    floor, ceiling = lh_bounds(length_heteroplasmy_threshold)
    if var_level < floor:
        return None, "BELOW_LH_FLOOR"
    naive_deleted = ref[len(var):]
    # `ref` is the VCF's own REF field - always the TRUE, un-baked genome,
    # regardless of what this sample's own major substitutions did.
    # `reference` is the (possibly baked - see main()'s pre-pass) genome
    # shift_deletion_right actually searches for repeats in. If a major/
    # near-homoplasmic substitution landed somewhere inside the deleted
    # span (e.g. T16189C at ~99%), the two disagree - the search would be
    # looking for a base (e.g. "T") that no longer exists at that
    # position in the reference it's searching, corrupting the result
    # (verified empirically on real data: HG01799's actual 16188 "CT->C"
    # call, with T16189C baked major, comes out as a nonsensical
    # "T16188t" using the naive REF-derived segment, versus the correct
    # "C16193c" using the baked one). Detected generally by comparing the
    # two directly, not by hardcoding the one position (16188/16189) this
    # was first caught at - the same staleness can happen at any baked
    # position (2026-09-23, user: "check whether 16188 hardcode has any
    # effect [...] could have for when true_mutect2 is used").
    del_start = pos + len(var) - 1
    baked_deleted = "".join(reference[del_start:del_start + len(naive_deleted)])
    deleted_segment = baked_deleted if baked_deleted != naive_deleted else naive_deleted
    pos, segment = shift_deletion_right(reference, pos, deleted_segment)

    is_major = var_level >= ceiling
    variant_parts = []
    for i, base in enumerate(segment):
        position = wrap_circular_position(pos + i)
        if is_major:
            variant_parts.append(f"{base}{position}-")
        else:
            variant_parts.append(f"{base}{position}{base.lower()}")

    updated_type = "DEL" if is_major else "LHP"
    return " ".join(variant_parts), updated_type

def format_variant(row, reference, min_variant_frequency, length_heteroplasmy_threshold):
    pos = int(row["Pos"])
    ref = row["Ref"]
    var = row["Variant"]
    var_type = row["Type"]
    var_level = row["VariantLevel"]

    if var_type == "SNP":
        return apply_snp(pos, ref, var, var_level, reference, min_variant_frequency)
    elif var_type == "INDEL":
        if len(ref) < len(var):
            return apply_insertion(pos, ref, var, var_level, reference, length_heteroplasmy_threshold)
        elif len(ref) > len(var):
            return apply_deletion(pos, ref, var, var_level, reference, length_heteroplasmy_threshold)

    return "N/A", var_type

def extract_numeric_value(empop_variant):
    match = re.search(r'\d+', empop_variant)
    return int(match.group()) if match else float('inf')

def extract_float_position(variant):
    match = re.search(r"(-?\d+\.?\d*)", str(variant))
    return float(match.group(1)) if match else float('inf')

def finalize_output_table(df, length_heteroplasmy_threshold):
    df["MUTECT2"] = df["MUTECT2"].astype(str).str.split()
    df = df.explode("MUTECT2").reset_index(drop=True)
    df["VariantLevel"] = pd.to_numeric(df["VariantLevel"], errors="coerce")

    def add_comma_separated_numbers(series):
        split_lists = series.dropna().astype(str).apply(lambda x: list(map(float, x.split(','))))
        if split_lists.empty:
            return ""
        summed = [sum(x) for x in zip(*split_lists)]
        return ",".join(f"{s:.4g}" for s in summed)

    # Rows are collapsed/summed when they carry the IDENTICAL exploded
    # label (e.g. "A16183M", or "-16193.1c") - the label already encodes
    # the full claim (position + resulting base/IUPAC code), so two rows
    # sharing it are always evidence for the exact same claim, regardless
    # of which original VCF record(s) they came from. That covers both:
    # different alleles of one multiallelic indel (a 1bp vs 2bp insertion,
    # both landing on "-16193.1c" after shifting), AND two independent
    # substitution events that both happen to assert "C at this position"
    # (e.g. a 2-base AA->CC substitution's second token and a separate
    # A->C call at the same position - both are "at least C here"
    # evidence and should combine). Rows with a DIFFERENT label - even at
    # the same numeric position, e.g. an unrelated indel's shifted label
    # landing near a SNP - are never summed, since they're different
    # claims; each stays its own row.
    def first_or_joined(series):
        # When a group combines multiple distinct source VCF records
        # (e.g. two alleles whose labels coincide, like "AA->CC"'s second
        # token and a separate "A->C" both producing "A16183M"), show all
        # of them rather than silently keeping only the first and
        # discarding the other's Pos/Ref/Variant/Filter - the combined
        # VariantLevel is correct, but its provenance shouldn't be hidden.
        vals = list(dict.fromkeys(series.dropna().astype(str)))
        if not vals:
            return None
        return vals[0] if len(vals) == 1 else " | ".join(vals)

    group_keys = ["MUTECT2"]
    numeric_agg = {
        "VariantLevel": "sum",
        "Coverage": add_comma_separated_numbers,
        "MeanBaseQuality": "first"
    }
    other_cols = [col for col in df.columns if col not in numeric_agg and col not in group_keys]
    full_agg = {**numeric_agg, **{col: first_or_joined for col in other_cols}}
    grouped = df.groupby(group_keys, as_index=False, dropna=False).agg(full_agg)

    grouped["variant_float_pos"] = grouped["MUTECT2"].apply(extract_float_position)
    grouped = grouped.sort_values("variant_float_pos").drop(columns=["variant_float_pos"])

    _, lh_ceiling = lh_bounds(length_heteroplasmy_threshold)

    def correct_length_het_case(row):
        if "." in row["MUTECT2"] and row["VariantLevel"] >= lh_ceiling:
            return row["MUTECT2"][:-1] + row["MUTECT2"][-1].upper()
        return row["MUTECT2"]

    grouped["MUTECT2"] = grouped.apply(correct_length_het_case, axis=1)
    return grouped

# bam_override mode (direct BAM read-counting substituted for Mutect2's
# own rows inside the boundary-run region, via call_region/call_
# boundary_run in call_repeat_regions.py) was removed here 2026-09-23
# (user: "it's not clean to overwrite mutect2 result with some new
# approach just to fit the output of fdstools. I think it's better to
# disable it if it provides junk and just rely on fdstools") - still
# available at git history if ever needed again. call_region and call_
# boundary_run themselves are untouched, still usable as a standalone
# tool via call_repeat_regions.py's own CLI.
def disable_homopolymer_length_calls(df, reference):
    """Drop Mutect2's own indel/length-axis rows inside MUTECT2_TARGET_
    REGIONS, with no replacement (the only remaining non-default mode,
    now that bam_override - which substituted direct BAM read counts -
    has been removed; see the comment above this function). Mutect2's
    ordinary substitution calls are left untouched, since an ordinary
    substitution never represents a length claim (see the same principle
    in report_separate_frame's docstring, process_fdstools_output_
    improved_better.py).

    Scoped to MUTECT2_TARGET_REGIONS - chrM:16180-16193 and chrM:303-315,
    matching the two regions process_fdstools_output_improved_better.py's
    own target_regions covers - this used to run genome-wide across all 9
    regions find_boundary_run_regions finds, on the reasoning that
    dropping a call needs no per-read reconciliation the way BAM-override
    does, so validation elsewhere didn't seem to matter. But the FDSTOOLS
    side was narrowed this same session, and leaving this side genome-
    wide meant a completely unrelated region (e.g. chrM:6416-6419,
    "A6419M") would show DISABLED on the Mutect2 side while FDSTOOLS
    reported it untouched - a real asymmetry, not two sides of one
    deliberate scope (2026-09-22, user: "I would restrict mutect2 as well
    same ways we do fdstools"; 2026-09-23, extended to 303-315 the same
    way: "I think we should extend it to 303-315 for fdstools and
    mutect2").

    Rows dropped here get no replacement row, and downstream (merge_
    fdstools_mutect2_improved.py, when homopolymer_mutect2_reporting is
    "disabled") their called_by_MUTECT2 is marked "DISABLED" rather than
    False, so a reader can't mistake "we didn't ask" for "Mutect2 asked
    and found nothing".
    """
    regions = MUTECT2_TARGET_REGIONS
    if not regions:
        return df

    df = df.copy()
    pos_numeric = pd.to_numeric(df["Pos"], errors="coerce")
    label_pos = df["MUTECT2"].apply(extract_float_position)
    is_length_type = df["Type"].astype(str).str.contains(r'\b(?:INS|DEL|LHP)\b', regex=True)

    drop_mask = pd.Series(False, index=df.index)
    for region in regions:
        leading = region["leading"]
        extension_start, extension_end = region["extension"]
        # A plain region (leading is None, e.g. 303-315) has no separate
        # reference-presence axis to drop - only its length-axis rows,
        # same as any position inside a true boundary-run's extension run.
        if leading:
            leading_start, leading_end = leading
            in_leading = pos_numeric.between(leading_start, leading_end) | label_pos.between(leading_start, leading_end)
        else:
            in_leading = pd.Series(False, index=df.index)
        # Not a plain between(extension_start, extension_end): a decimal
        # insertion label anchored at the extension run's own far end
        # (e.g. "-16193.1c") extracts to 16193.1, just past extension_end -
        # between() would miss it entirely.
        in_extension = (
            pos_numeric.between(extension_start, extension_end)
            | ((label_pos >= extension_start) & (label_pos < extension_end + 1))
        )
        drop_mask |= in_leading | (is_length_type & in_extension)

    for _, dropped in df.loc[drop_mask].iterrows():
        print(f"NOTE: dropping Mutect2 row {dropped['MUTECT2']!r} (VariantLevel="
              f"{dropped['VariantLevel']}) - homopolymer_mutect2_reporting=disabled, "
              f"no replacement", file=sys.stderr)

    return df.loc[~drop_mask]


def main():
    parser = argparse.ArgumentParser(description="Process mitochondrial variants into EMPOP format.")
    parser.add_argument("input_file", help="Input TSV file with variants")
    parser.add_argument("output_file", help="Output TSV file")
    parser.add_argument("reference_fasta", help="Reference genome in FASTA format")
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minor allele frequency threshold")
    parser.add_argument("--lh_thresh", type=float, default=10.0, help="Symmetric length-heteroplasmy threshold (floor and, via 1-threshold, ceiling), e.g. 10.0 -> report only 10-90%%, lowercase in between, major above 90%%")
    parser.add_argument("--homopolymer_mutect2_reporting", choices=["disabled", "true_mutect2"], default="true_mutect2",
                         help="How Mutect2's own calls are handled inside the boundary-run region (leading run + different-base extension run, chrM:16180-16193 - "
                              "the only one either mode touches; every other boundary-run region genome-wide is left as Mutect2 reported it). "
                              "disabled drops Mutect2's indel/length-axis rows there with no replacement, keeping ordinary substitution calls like T16189C untouched (see "
                              "disable_homopolymer_length_calls); merge_fdstools_mutect2_improved.py marks called_by_MUTECT2 as DISABLED for these when this mode is active. "
                              "true_mutect2 (default) leaves Mutect2's own calls as is. (A third mode, bam_override, substituted direct BAM read-counts for Mutect2's "
                              "rows there - removed 2026-09-23 as not clean to overwrite Mutect2's own results just to fit FDSTOOLS' output shape; use disabled and "
                              "rely on FDSTOOLS for this region instead.)")

    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_file, sep="\t")
        reference = load_reference_list(args.reference_fasta)
        # Unbaked snapshot, passed to disable_homopolymer_length_calls
        # below - kept distinct from the baked `reference` used for the
        # indel repeat-shift searches right below, in case a future
        # homopolymer_mutect2_reporting mode needs the true, un-baked
        # reference the way the removed bam_override mode did.
        reference_unbaked = list(reference)

        # Bake MAJOR (near-homoplasmic) substitutions into the reference
        # used for indel repeat-shift searches, so a dominant nearby
        # substitution (e.g. T16189C at 98.7%) is correctly treated as
        # part of the observed sequence when shifting an insertion/
        # deletion through it - the C-stretch really does extend through
        # a near-homoplasmic T->C there. A minority/heteroplasmic
        # substitution is deliberately NOT baked in, since most reads
        # still show the true reference at that position. This is a
        # single order-independent pre-pass (each major substitution only
        # touches its own position(s)), not an order-dependent side
        # effect of per-row processing.
        major_snp_threshold = 1 - (args.min_vf / 100)
        for _, snp_row in df.iterrows():
            if str(snp_row.get("Type")) != "SNP":
                continue
            level = pd.to_numeric(snp_row.get("VariantLevel"), errors="coerce")
            if pd.isna(level) or level < major_snp_threshold:
                continue
            ref_str, var_str = str(snp_row["Ref"]), str(snp_row["Variant"])
            if len(ref_str) != len(var_str):
                continue
            snp_pos = int(snp_row["Pos"])
            for i, v in enumerate(var_str):
                reference[snp_pos + i - 1] = v

        results = []
        types = []
        for idx in reversed(df.index):
            row = df.loc[idx]
            variant, updated_type = format_variant(row, reference, args.min_vf/100, args.lh_thresh/100)
            results.append(variant)
            types.append(updated_type)

        df["MUTECT2"] = results[::-1]
        df["Type"] = types[::-1]

        # Drop length variants below the symmetric length-heteroplasmy floor entirely
        df = df[df["Type"] != "BELOW_LH_FLOOR"].reset_index(drop=True)

        df = finalize_output_table(df, args.lh_thresh/100)

        if args.homopolymer_mutect2_reporting == "disabled":
            df = disable_homopolymer_length_calls(df, reference_unbaked)

        # Rename selected columns
        df = df.rename(columns={
            "VariantLevel": "vf_MT2",
            "Coverage": "rd_MT2",
            "MeanBaseQuality": "MBQ"
        })

        # Move 'MUTECT2' column to the front
        cols = df.columns.tolist()
        if "MUTECT2" in cols:
            cols.insert(0, cols.pop(cols.index("MUTECT2")))
            df = df[cols]
            
        df.to_csv(args.output_file, sep="\t", index=False)
        print(f"Processed file saved to {args.output_file}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
