#!/usr/bin/env python3

import pandas as pd
import re
import argparse
import sys
import os
import traceback

from Bio import SeqIO

import frames
import iupac
import notation

# Utility: extract numeric position from sequence
def extract_position(seq):
    match = re.search(r"(\d+\.?\d*)", seq)
    return float(match.group(1)) if match else float('inf')

# Adjust for circular mtDNA positions (16570–16587 → 1–18)
def adjust_circular_position(pos):
    if 16570 <= pos <= 16587:
        return pos - 16569
    return pos

LH_DROP_SENTINEL = "__BELOW_LH_FLOOR__"

def lh_bounds_pct(threshold_pct):
    """Symmetric floor/ceiling (percentage scale, 0-100) around a single
    length-heteroplasmy threshold, e.g. threshold_pct=10 -> (10, 90).
    Accepts either side (10 or 90) and always returns (floor, ceiling)
    with floor <= ceiling."""
    return min(threshold_pct, 100 - threshold_pct), max(threshold_pct, 100 - threshold_pct)

# Length heteroplasmy: minor insertions and deletions in lowercase
def resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold):
    seq = row['sequence']

    if 'DEL' in seq or '.' in seq:
        lh_floor, lh_ceiling = lh_bounds_pct(length_heteroplasmy_threshold)
        if row['variant_frequency'] < lh_floor:
            return LH_DROP_SENTINEL
        is_major = row['variant_frequency'] >= lh_ceiling
        # Minor deletion as ref base in lowercase (A523a), matching the Mutect2 side
        if 'DEL' in seq:
            return seq.replace('DEL', '-' if is_major else seq[0].lower())
        return '-' + seq if is_major else '-' + seq[:-1] + seq[-1].lower()
    return seq


SUBSTITUTION = re.compile(r"^[ACGT]\d+[ACGT]$")


def position_rows(final, other_share, min_vf):
    """One row per position for FDSTOOLS' substitution labels (C756A, C756T), by
    iupac.call: the bases with at least min_vf are present, rCRS when what
    other_share[position] (every substitution and deletion there, also those under
    min_vf) leaves reaches min_vf. Several bases other than rCRS give one row with
    frequency and reads per base ("A 30, T 10")."""
    labels, rows, combined = {}, [], []
    subs = final[final["sequence"].str.match(SUBSTITUTION)]
    for pos, group in subs.groupby("position"):
        ref = group["sequence"].iloc[0][0]
        shares = dict(zip(group["sequence"].str[-1], group["variant_frequency"]))
        code, others = iupac.call(ref, {**shares, ref: 100 - other_share.get(pos, 0)}, min_vf)
        label = f"{ref}{int(pos)}{code}"
        if len(group) == 1:
            labels[group.index[0]] = label
            continue
        reads = dict(zip(group["sequence"].str[-1], group["total"]))
        first = group.loc[group["variant_frequency"].idxmax()].to_dict()
        rows.append({**first, "sequence": label, "variant_frequency": iupac.format_values(others),
                     "total": iupac.format_values({b: reads[b] for b in others}, 0),
                     "interpolated_total_coverage": group["interpolated_total_coverage"].max(),
                     "num_markers": group["num_markers"].max()})
        combined += list(group.index)
    final = final.copy()
    for i, label in labels.items():
        final.at[i, "sequence"] = label
    if rows:
        final = pd.concat([final.drop(index=combined), pd.DataFrame(rows)], ignore_index=True)
    return final

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

def in_window(label, region):
    """Whether an FDSTOOLS label (A16183C, T16189DEL, 16193.1C) lies in a frame region or
    within frames.FLANK bases of it, where the rows come from the frames' alignment."""
    m = re.match(r"^[ACGTN]?(\d+)(?:\.\d+)?([ACGT]|DEL)?$", label)
    return bool(m) and region.first - frames.FLANK <= int(m.group(1)) <= region.last + frames.FLANK


def frame_rows(tssv_path, reference, marker_map, marker_total_reads, frame, min_vf, lh_thresh):
    """Rows for 57-60, 300-315 and 16180-16193 from the frames: every sequence of the
    amplicon holding a region (tssv.csv) is placed and written in the chosen frame, out
    of the same reads as every other row of that amplicon. The changes within
    frames.FLANK bases of the region come from the same alignment, so a change next to
    the region is never counted on both sides of its edge. The region's rows also get
    the dominant molecule's change at their position (frames.dominant_cell)."""
    tssv = pd.read_csv(tssv_path, sep="\t", dtype=str)
    ranges = {m: tuple(int(x) for x in r.split("-")) for m, r in marker_map.items()}
    rows = []
    for region in frames.REGIONS.values():
        for marker in frames.covering(region, ranges):
            coverage = marker_total_reads.get(marker, 0)
            if not coverage:
                continue
            start, end = ranges[marker]
            amplicon = tssv[tssv["marker"] == marker]
            sequences = list(zip(amplicon["sequence"], pd.to_numeric(amplicon["total"])))
            molecules = frames.region_molecules(sequences, start, end, region, reference)
            region_rows = frames.rows(region, molecules, coverage, frame, reference, min_vf, lh_thresh)
            flank = frames.flank_rows(sequences, start, end, region, reference, coverage, min_vf, lh_thresh)
            dominant = frames.dominant(region, molecules, frame, reference)
            for label, share, note in [r + (f"{frame} frame",) for r in region_rows] + [r + (None,) for r in flank]:
                cell = frames.dominant_cell(label, dominant, coverage)
                if isinstance(share, dict):  # several bases other than rCRS
                    total = iupac.format_values({b: s * coverage / 100 for b, s in share.items()}, 0)
                    share = iupac.format_values(share)
                else:
                    total, share = round(share * coverage / 100), round(share, 2)
                rows.append({"sequence": label, "total": total,
                             "interpolated_total_coverage": coverage, "variant_frequency": share,
                             "marker": marker, "num_markers": 1, "is_noise_or_low_frq": False,
                             "variant_note": note, "dominant_molecule": cell or None})
    return pd.DataFrame(rows)


def replace_labels(table, old, new, note):
    """Rows old {index: label} rewritten as the labels new: a row whose label is among the
    new ones keeps its frequency, the other rows go, and each new label left over gets
    the lowest frequency of the rows that went, with `note` in variant_note."""
    staying = {i: label for i, label in old.items() if label in new}
    going = [i for i in old if i not in staying]
    added = [label for label in new if label not in staying.values()]
    rows = []
    if added:
        lowest = table.loc[going or list(old)].sort_values("variant_frequency").iloc[0].to_dict()
        rows = [{**lowest, "sequence": label, "variant_note": note} for label in added]
    return pd.concat([table.drop(index=going), pd.DataFrame(rows)], ignore_index=True)


def respell_rows(table, reference):
    """Groups of major labels outside the frame regions written again by the general
    rule (notation.respell_majors), so both callers spell a molecule alike."""
    protected = [(r.first, r.last) for r in frames.REGIONS.values()]
    for old, new in notation.respell_majors(table["sequence"], reference, protected):
        rows = {i: label for i, label in table["sequence"].items() if label in old}
        table = replace_labels(table, rows, new, f"written from {' '.join(old)}")
    return table


# Main processing function
def process_fdstools_sast(file_path, marker_map_path, output_file, min_variant_frequency_pct=5.0, depth_threshold=10, length_heteroplasmy_threshold=90.0,
                          tssv_path=None, reference_path=None, frame="separate"):
    df = pd.read_csv(file_path, sep="\t", dtype=str)
    df = df.drop(columns=[
        "total_mp_max", "forward_pct", "forward", "forward_mp_sum",
        "forward_mp_max", "reverse", "reverse_mp_sum", "reverse_mp_max"
    ], errors="ignore")

    df["total_mp_sum"] = pd.to_numeric(df["total_mp_sum"], errors="coerce").fillna(0)
    df["total"] = pd.to_numeric(df["total"], errors="coerce").fillna(0)

    # LOW: amplicon whose reads, summed over all its rows except "Other sequences"
    # (the same reads the frequencies are computed on), are below depth_threshold.
    # Checked for every amplicon in the library, so a complete dropout is flagged
    # even if it has no row in the file.
    marker_map = load_marker_ranges(marker_map_path)
    marker_total_reads = df[df["sequence"] != "Other sequences"].groupby("marker")["total"].sum()
    single_low_coverage = pd.DataFrame([
        {"marker": marker, "total": marker_total_reads.get(marker, 0), "sequence": "LOW"}
        for marker in marker_map
        if marker_total_reads.get(marker, 0) < depth_threshold
    ], columns=["marker", "total", "sequence"])


    df["total"] = df["total"].fillna(0)
    df["total_mp_sum"] = df["total_mp_sum"].fillna(0)
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

    grouped = df.groupby(["marker", "sequence"], as_index=False).agg(
        total=("total", "sum"),
        total_mp_sum=("total_mp_sum", "sum"),
    )

    # Coverage: reads counted per amplicon, excluding "Other sequences" (the same
    # count LOW uses); clipped to 1 to avoid division by zero
    grouped["interpolated_total_coverage"] = grouped["marker"].map(marker_total_reads).fillna(0).clip(lower=1)


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

    # In the frame regions the rows come from the frames instead of FDSTOOLS' labels
    rows_from_frames = pd.DataFrame()
    reference = str(SeqIO.read(reference_path, "fasta").seq) if reference_path else None
    if tssv_path:
        rows_from_frames = frame_rows(tssv_path, reference, marker_map, marker_total_reads, frame,
                                      min_variant_frequency_pct, length_heteroplasmy_threshold)
        in_a_window = final["sequence"].map(
            lambda label: any(in_window(label, region) for region in frames.REGIONS.values())).astype(bool)
        final = final[~in_a_window]

    drop_seqs = ["Other", "sequences", "REF", "N3107DEL", "No", "data"]
    final = final[(~final["sequence"].isin(drop_seqs))]

    # Share of the molecules with another base or a deletion at each position, all of
    # them (also under min_vf): what is left is rCRS
    changed = final["sequence"].str.match(r"^[ACGT]\d+(?:[ACGT]|DEL)$")
    other_share = final[changed].groupby("position")["variant_frequency"].sum().to_dict()

    final["is_noise_or_low_frq"] = (final["sequence"].isin(["Other sequences"])) | (final["variant_frequency"] < min_variant_frequency_pct)
    final = final[~final["is_noise_or_low_frq"]]

    # If nothing remains after filtering (e.g., H2O / No data), write empty output and stop
    if final.empty and rows_from_frames.empty:
        pd.DataFrame(columns=[
            "FDSTOOLS", "vf_FDS", "rd_FDS", "interpolated_total_coverage",
            "variant_note", "marker", "marker_range", "num_markers"
        ]).to_csv(output_file, sep="\t", index=False)
        print(f"Output written to {output_file} (no variants)")
        return
    if not final.empty:
        final = position_rows(final, other_share, min_variant_frequency_pct)
        resolved = final.apply(
            lambda row: resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold),
            axis=1
        )

        # Force to a plain 1D Series of strings
        final["sequence"] = pd.Series(resolved, index=final.index).astype(str)

        # Drop length variants below the symmetric length-heteroplasmy floor entirely
        final = final[final["sequence"] != LH_DROP_SENTINEL]

        if reference is not None:
            final = respell_rows(final, reference)
    
    # final["sequence"] = final.apply(
    #     lambda row: resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold),
    #     axis=1
    # )

    final = pd.concat([final, rows_from_frames, single_low_coverage], ignore_index=True)

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
        "dominant_molecule",
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
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minimum variant frequency threshold")
    parser.add_argument("--depth", type=int, default=10, help="Read depth threshold for low coverage")
    parser.add_argument("--lh_thresh", type=float, default=10.0, help="Symmetric length-heteroplasmy threshold (floor and, via 100-threshold, ceiling), e.g. 10.0 -> report only 10-90%%, lowercase in between, major above 90%%")
    parser.add_argument("--tssv", help="FDSTOOLS tssv.csv (sequences per amplicon); with it, 57-60, 300-315 and 16180-16193 are written in the frames")
    parser.add_argument("--reference", help="rCRS_NimaGen.fasta; with it, major calls outside the frame regions are spelled by the general rule (needed with --tssv)")
    parser.add_argument("--frame", choices=["separate", "shared"], default="separate", help="Frame for 16180-16193 and 300-315")
    args = parser.parse_args()

    try:
        process_fdstools_sast(
            file_path=args.input,
            marker_map_path=args.marker_map,
            output_file=args.output,
            min_variant_frequency_pct=args.min_vf,
            depth_threshold=args.depth,
            length_heteroplasmy_threshold=args.lh_thresh,
            tssv_path=args.tssv,
            reference_path=args.reference,
            frame=args.frame
        )
    except Exception:
        traceback.print_exc()
        sys.exit(1)
    # except Exception as e:
    #     print(f"Error: {e}", file=sys.stderr)
    #     sys.exit(1)

if __name__ == "__main__":
    main()
