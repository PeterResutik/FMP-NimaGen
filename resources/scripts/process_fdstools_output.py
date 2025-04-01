import argparse
import pandas as pd # type: ignore
import re

# Step 10: Extract numeric positions from the sequence column for sorting
def extract_position(seq):
    match = re.search(r"(\d+\.?\d*)", seq)
    return float(match.group(1)) if match else float('inf')

# Define function to process FDSTools SAST TSV files
def process_fdstools_sast(file_path: str, threshold: float = 8.0):
    IUPAC_CODES = {
        frozenset(["A", "G"]): "R",
        frozenset(["C", "T"]): "Y",
        frozenset(["A", "C"]): "M",
        frozenset(["G", "T"]): "K",
        frozenset(["G", "C"]): "S",
        frozenset(["A", "T"]): "W"
    }

    # Load the FDSTools SAST TSV file
    df = pd.read_csv(file_path, sep="\t", dtype=str)
    df["total_mp_sum"] = pd.to_numeric(df["total_mp_sum"], errors="coerce")
    df["total"] = pd.to_numeric(df["total"], errors="coerce")

    # Step 1: Calculate total read depth per marker
    df["total"] = df["total"].fillna(0)
    df["total_mp_sum"] = df["total_mp_sum"].fillna(0)

    # Step 2: Recalculate total read depth per marker excluding noise and low-confidence
    df["is_noise_or_low"] = (df["sequence"].isin(["Other sequences"])) | (df["total_mp_sum"] < threshold)

    columns_to_drop = [
    "total_mp_max", "forward_pct", "forward", "forward_mp_sum",
    "forward_mp_max", "reverse", "reverse_mp_sum", "reverse_mp_max"
    ]

    df = df.drop(columns=columns_to_drop, errors="ignore")
    
    clean_total_per_marker = df[~df["is_noise_or_low"]].groupby("marker")["total"].sum().rename("clean_marker_total_wo_OS_THR")
    df = df.merge(clean_total_per_marker, on="marker", how="left")

    # Step 3: Compute normalized variant frequency (only for retained rows)
    df["variant_frequency_wo_OS_THR"] = (df["total"] / df["clean_marker_total_wo_OS_THR"] * 100).round(2)

    # Step 5: Split multiple variants
    df = df.assign(sequence=df["sequence"].str.split())
    df = df.explode("sequence").reset_index(drop=True)
    
    # Step 4: Flag rows to retain
    drop_seqs = ["Other", "sequences", "REF", "N3107DEL"]
    df = df[(~df["sequence"].isin(drop_seqs)) & (df["total_mp_sum"] >= threshold)].copy()


    # Step 6: Calculate estimated coverage for each row
    df["estimated_total_coverage"] = (
        df["total"] / (df["total_mp_sum"] / 100)
    ).round(0).astype("Int64")


    # Step 7: Group by marker + sequence to sum within same marker
    grouped_same_marker = df.groupby(["marker", "sequence"], as_index=False).agg(
        # total=("total", "sum"),
        # total_mp_sum=("total_mp_sum", "sum"),
        # estimated_total_coverage=("estimated_total_coverage", "sum")
        total=("total", "sum"),
        total_mp_sum=("total_mp_sum", "sum"),
        estimated_total_coverage=("estimated_total_coverage", "max"),
        is_noise_or_low=("is_noise_or_low", "first"),
        clean_marker_total_wo_OS_THR=("clean_marker_total_wo_OS_THR", "first"),
        variant_frequency_wo_OS_THR=("variant_frequency_wo_OS_THR", "sum")
    )
    # # print(grouped_same_marker)
    # # Step 8: Recompute variant_frequency within marker
    # grouped_same_marker["variant_frequency"] = (
    #     grouped_same_marker["total"] / grouped_same_marker["estimated_total_coverage"] * 100
    # ).round(1)

    # print(grouped_same_marker)
    
    # Step 9: Group across markers to merge overlapping amplicons (same variant)
    # df["estimated_total_coverage_across_markers"] = (df["total"] / (df["total_mp_sum"] / 100)).round(0).astype("Int64")
    grouped_final = grouped_same_marker.groupby("sequence", as_index=False).agg(
        total=("total", "sum"),
        # total_mp_sum=("total_mp_sum", "sum"),
        estimated_total_coverage=("estimated_total_coverage", "sum"),
        # estimated_total_coverage_across_markers=("estimated_total_coverage_across_markers", "sum"),
        is_noise_or_low=("is_noise_or_low", "first"),
        clean_marker_total_wo_OS_THR=("clean_marker_total_wo_OS_THR", "sum"),
        # variant_frequency_wo_OS_THR=("variant_frequency_wo_OS_THR", "sum"),
        num_markers=("marker", "nunique")
    )


    # print(grouped_final)
    
    grouped_final["variant_frequency"] = (
        grouped_final["total"] / grouped_final["estimated_total_coverage"] * 100
    ).round(1)

    grouped_final["variant_frequency_wo_OS_THR"] = (
        grouped_final["total"] / grouped_final["clean_marker_total_wo_OS_THR"] * 100
    ).round(1)

    grouped_final["position"] = grouped_final["sequence"].apply(extract_position)
    grouped_final = grouped_final.sort_values(by="position").drop(columns=["position"])

    
    # Step 11: Apply IUPAC codes for heteroplasmies
    # Step 11: Apply IUPAC codes and adjust formatting for heteroplasmies
    def resolve_heteroplasmy(row):
        seq = row['sequence']
    
        # Handle deletions
        if 'DEL' in seq:
            return seq.replace('DEL', '-')

        # Handle insertions: prefix with "-"

        # Handle insertions / length heteroplasmies
        if '.' in seq:
            if row['variant_frequency_wo_OS_THR'] < 92:
                return '-' + seq[:-1] + seq[-1].lower()  # e.g. 309.2C -> 309.2c
            else:
                return '-' + seq # leave as-is if frequency is high

        # Handle point heteroplasmies with IUPAC
        if row['variant_frequency_wo_OS_THR'] < 92:
            match = re.match(r'([ACGT])(\d+)([ACGT])', seq)
            if not match:
                return seq
            ref, pos, alt = match.groups()
            code = IUPAC_CODES.get(frozenset([ref, alt]))
            if code:
                return f"{ref}{pos}{code}"
    
        return seq



    grouped_final['sequence'] = grouped_final.apply(resolve_heteroplasmy, axis=1)

    # print(grouped_final)
    
    return grouped_final
    # return grouped_final, df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process FDSTools SAST variant table.")
    parser.add_argument("input_file", help="Path to the input SAST TSV file.")
    parser.add_argument("output_file", help="Path to save the processed output TSV.")
    parser.add_argument("--threshold", type=float, default=8.0, help="Threshold for filtering total_mp_sum (default: 8.0)")

    args = parser.parse_args()

    final_df = process_fdstools_sast(args.input_file, threshold=args.threshold)
    final_df.to_csv(args.output_file, sep="\t", index=False)
    # print(f"Processed variant table saved to: {args.output_file}")

