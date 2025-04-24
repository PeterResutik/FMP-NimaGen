#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import argparse
import sys
import os

# Utility: extract numeric position from sequence
def extract_position(seq):
    match = re.search(r"(\d+\.?\d*)", seq)
    return float(match.group(1)) if match else float('inf')

# Adjust for circular mtDNA positions (16570–16587 → 1–18)
def adjust_circular_position(pos):
    if 16570 <= pos <= 16587:
        return pos - 16569
    return pos

# IUPAC resolution for heteroplasmies
def resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold, IUPAC_CODES):
    seq = row['sequence']

    if 'DEL' in seq:
        return seq.replace('DEL', '-')
    if '.' in seq:
        if row['variant_frequency_wo_noise_or_low_frq'] < length_heteroplasmy_threshold:
            return '-' + seq[:-1] + seq[-1].lower()
        else:
            return '-' + seq
    if row['variant_frequency_wo_noise_or_low_frq'] < 100 - min_variant_frequency_pct:
        match = re.match(r'([ACGT])(\d+)([ACGT])', seq)
        if match:
            ref, pos, alt = match.groups()
            code = IUPAC_CODES.get(frozenset([ref, alt]))
            if code:
                return f"{ref}{pos}{code}"
    return seq

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
def process_fdstools_sast(file_path, marker_map_path, output_file, min_variant_frequency_pct=5.0, depth_threshold=10, length_heteroplasmy_threshold=90.0):
    IUPAC_CODES = {
        frozenset(["A", "G"]): "R", frozenset(["C", "T"]): "Y",
        frozenset(["A", "C"]): "M", frozenset(["G", "T"]): "K",
        frozenset(["G", "C"]): "S", frozenset(["A", "T"]): "W"
    }

    df = pd.read_csv(file_path, sep="\t", dtype=str)
    df = df.drop(columns=[
        "total_mp_max", "forward_pct", "forward", "forward_mp_sum",
        "forward_mp_max", "reverse", "reverse_mp_sum", "reverse_mp_max"
    ], errors="ignore")

    df["total_mp_sum"] = pd.to_numeric(df["total_mp_sum"], errors="coerce").fillna(0)
    df["total"] = pd.to_numeric(df["total"], errors="coerce").fillna(0)

    marker_counts = df["marker"].value_counts()
    single_low_coverage = df[
        df["marker"].isin(marker_counts[marker_counts == 1].index) & 
        (df["total"] < depth_threshold)
    ].copy()
    single_low_coverage["sequence"] = "LOW OR NO COVERAGE"
    single_low_coverage = single_low_coverage.drop(columns=["total_mp_sum", "flags"], errors="ignore")

    df["is_noise_or_low_frq"] = df["sequence"].isin(["Other sequences"]) | (df["total_mp_sum"] < min_variant_frequency_pct)
    clean_total_per_marker = df[~df["is_noise_or_low_frq"]].groupby("marker")["total"].sum().rename("total_wo_noise_or_low_frq")
    df = df.merge(clean_total_per_marker, on="marker", how="left")
    df["variant_frequency_wo_noise_or_low_frq"] = (df["total"] / df["total_wo_noise_or_low_frq"] * 100).round(2)

    df = df.assign(sequence=df["sequence"].str.split()).explode("sequence").reset_index(drop=True)
    drop_seqs = ["Other", "sequences", "REF", "N3107DEL"]
    df = df[(~df["sequence"].isin(drop_seqs)) & (df["total_mp_sum"] >= min_variant_frequency_pct)].copy()
    df["interpolated_total_coverage"] = (np.ceil(df["total"] / (df["total_mp_sum"] / 100))).astype("Int64")

    grouped = df.groupby(["marker", "sequence"], as_index=False).agg(
        total=("total", "sum"),
        total_mp_sum=("total_mp_sum", "sum"),
        interpolated_total_coverage=("interpolated_total_coverage", "max"),
        is_noise_or_low_frq=("is_noise_or_low_frq", "first"),
        total_wo_noise_or_low_frq=("total_wo_noise_or_low_frq", "first"),
        variant_frequency_wo_noise_or_low_frq=("variant_frequency_wo_noise_or_low_frq", "sum")
    )

    final = grouped.groupby("sequence", as_index=False).agg(
        marker=("marker", "first"),
        total=("total", "sum"),
        interpolated_total_coverage=("interpolated_total_coverage", "sum"),
        is_noise_or_low_frq=("is_noise_or_low_frq", "first"),
        total_wo_noise_or_low_frq=("total_wo_noise_or_low_frq", "sum"),
        num_markers=("marker", "nunique")
    )

    final["variant_frequency"] = (final["total"] / final["interpolated_total_coverage"] * 100).round(2)
    final["variant_frequency_wo_noise_or_low_frq"] = (final["total"] / final["total_wo_noise_or_low_frq"] * 100).round(2)
    final["position"] = final["sequence"].apply(extract_position)

    # Merge substitutions and deletions at same position
    merged_rows, used_indices = [], set()
    for pos, group in final.groupby("position"):
        if group.shape[0] != 2: continue
        del_row = group[group["sequence"].str.endswith("DEL")]
        sub_row = group[~group["sequence"].str.endswith("DEL")]
        if del_row.empty or sub_row.empty: continue
        del_idx, sub_idx = del_row.index[0], sub_row.index[0]
        if del_idx in used_indices or sub_idx in used_indices: continue

        total = del_row["total"].iloc[0] + sub_row["total"].iloc[0]
        coverage = del_row["total_wo_noise_or_low_frq"].iloc[0]
        freq = round(total / coverage * 100, 1) if coverage else 0
        ref, pos_str, alt = re.match(r'([ACGT])(\d+)([ACGT])', sub_row["sequence"].iloc[0]).groups()
        del_freq = del_row["variant_frequency_wo_noise_or_low_frq"].iloc[0]
        merged_seq = f"{ref}{pos_str}{alt.lower()}" if del_freq < length_heteroplasmy_threshold else f"{ref}{pos_str}{alt}"
        freq_annotation = (
            f"sub:{sub_row['variant_frequency'].iloc[0]} ({sub_row['variant_frequency_wo_noise_or_low_frq'].iloc[0]}) | "
            f"del:{del_row['variant_frequency'].iloc[0]} ({del_row['variant_frequency_wo_noise_or_low_frq'].iloc[0]})"
        )

        merged_rows.append({
            "sequence": merged_seq,
            "total": total,
            "interpolated_total_coverage": del_row["interpolated_total_coverage"].iloc[0],
            "is_noise_or_low_frq": False,
            "total_wo_noise_or_low_frq": total,
            "num_markers": f"sub:{sub_row['num_markers'].iloc[0]} | del:{del_row['num_markers'].iloc[0]}",
            "variant_frequency": freq_annotation,
            "variant_frequency_wo_noise_or_low_frq": freq,
            "marker": sub_row["marker"].iloc[0],
            "position": float(pos)
        })

        used_indices.update([del_idx, sub_idx])

    final = final.drop(index=used_indices)
    if merged_rows:
        final = pd.concat([final, pd.DataFrame(merged_rows)], ignore_index=True)

    final["sequence"] = final.apply(
        lambda row: resolve_heteroplasmy(row, min_variant_frequency_pct, length_heteroplasmy_threshold, IUPAC_CODES),
        axis=1
    )

    final = pd.concat([final, single_low_coverage], ignore_index=False)

    marker_map = load_marker_ranges(marker_map_path)
    final["marker_range"] = final["marker"].map(marker_map)
    final["position"] = final["marker"].apply(extract_position)
    final["position"] = final["position"].apply(adjust_circular_position)
    final = final.sort_values(by="position").drop(columns=["position"])

    final.to_csv(output_file, sep="\t", index=False)
    print(f"✔ Output written to {output_file}")

# CLI wrapper
def main():
    parser = argparse.ArgumentParser(description="Process FDSTools SAST TSV file with heteroplasmy handling.")
    parser.add_argument("input", help="Input SAST TSV file")
    parser.add_argument("output", help="Output processed file")
    parser.add_argument("--marker_map", help="Path to marker map file")
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minimum variant frequency threshold")
    parser.add_argument("--depth", type=int, default=10, help="Read depth threshold for low coverage")
    parser.add_argument("--lh_thresh", type=float, default=90.0, help="Length heteroplasmy frequency threshold")
    args = parser.parse_args()

    try:
        process_fdstools_sast(
            file_path=args.input,
            marker_map_path=args.marker_map,
            output_file=args.output,
            min_variant_frequency_pct=args.min_vf,
            depth_threshold=args.depth,
            length_heteroplasmy_threshold=args.lh_thresh
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
