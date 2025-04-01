import argparse
import pandas as pd # type: ignore
import re

def process_empop_variant_table(file_path: str) -> pd.DataFrame:
    """
    Processes a variant table with EMPOP-style variant annotations.
    - Splits multi-variant rows
    - Sums VariantLevel and allele-specific Coverage
    - Keeps the first value of MeanBaseQuality
    - Sorts variants by position

    Parameters:
    - file_path: path to the tab-separated file

    Returns:
    - A cleaned and sorted DataFrame
    """

    # Load file
    df = pd.read_csv(file_path, sep="\t")

    # Split EMPOP_Variant column and explode into rows
    df["EMPOP_Variant"] = df["EMPOP_Variant"].astype(str).str.split()
    df = df.explode("EMPOP_Variant").reset_index(drop=True)

    # Convert numeric columns
    df["VariantLevel"] = pd.to_numeric(df["VariantLevel"], errors="coerce")

    # Helper: sum comma-separated numbers elementwise
    def add_comma_separated_numbers(series):
        split_lists = series.dropna().astype(str).apply(lambda x: list(map(float, x.split(','))))
        if split_lists.empty:
            return ""
        summed = [sum(x) for x in zip(*split_lists)]
        return ",".join(f"{s:.4g}" for s in summed)

    # Aggregate
    group_keys = ["EMPOP_Variant"]
    numeric_agg = {
        "VariantLevel": "sum",
        "Coverage": add_comma_separated_numbers,
        "MeanBaseQuality": "first"
    }
    other_cols = [col for col in df.columns if col not in numeric_agg and col not in group_keys]
    full_agg = {**numeric_agg, **{col: "first" for col in other_cols}}

    grouped = df.groupby(group_keys, as_index=False).agg(full_agg)

    # Sort by numeric position extracted from variant
    def extract_position(variant):
        match = re.search(r"(\d+\.?\d*)", str(variant))
        return float(match.group(1)) if match else float('inf')

    grouped["position"] = grouped["EMPOP_Variant"].apply(extract_position)
    grouped = grouped.sort_values(by="position").drop(columns=["position"])

    def correct_length_het_case(row):
        if "." in row["EMPOP_Variant"] and row["VariantLevel"] >= 0.92:
            return row["EMPOP_Variant"][:-1] + row["EMPOP_Variant"][-1].upper()
        return row["EMPOP_Variant"]

    grouped["EMPOP_Variant"] = grouped.apply(correct_length_het_case, axis=1)
    
    return grouped

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process EMPOP-style variant table.")
    parser.add_argument("input_file", help="Path to input variant file (TSV format)")
    parser.add_argument("output_file", help="Path to output processed file (TSV format)")
    args = parser.parse_args()

    df_processed = process_empop_variant_table(args.input_file)
    df_processed.to_csv(args.output_file, sep="\t", index=False)
    # print(f"Processed EMPOP variant table saved to: {args.output_file}")

