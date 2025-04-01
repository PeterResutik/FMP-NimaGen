import argparse
import pandas as pd # type: ignore
import re

def merge_variant_callers(file1: str, file2: str) -> pd.DataFrame:
    """
    Merges two variant tables from different callers on their variant column.
    
    Parameters:
        file1: Path to the first variant caller table (expects 'sequence' column)
        file2: Path to the second variant caller table (expects 'EMPOP_Variant' column)
        
    Returns:
        A merged DataFrame with flags and full column preservation, sorted by position.
    """
    # Load both files
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")
    
    # Rename variant columns to common key
    df1 = df1.rename(columns={"sequence": "variant"})
    df2 = df2.rename(columns={"EMPOP_Variant": "variant"})
    
    # Merge the dataframes on the variant column
    merged = pd.merge(df1, df2, on="variant", how="outer", suffixes=("_vc1", "_vc2"))

    # Add flags for presence in each caller
    merged["called_in_vc1"] = ~merged["total"].isna()
    merged["called_in_vc2"] = ~merged["VariantLevel"].isna()

    # Extract numeric position for sorting
    def extract_position(seq):
        match = re.search(r"(\d+\.?\d*)", str(seq))
        return float(match.group(1)) if match else float('inf')

    merged["variant_position"] = merged["variant"].apply(extract_position)
    merged = merged.sort_values(by="variant_position").drop(columns=["variant_position"])

    # Reorder columns for clarity
    front = ["variant", "called_in_vc1", "called_in_vc2"]
    other = [col for col in merged.columns if col not in front]
    merged = merged[front + other]

    return merged

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge mitochondrial variant calls from two sources.")
    parser.add_argument("caller1", help="Path to the first variant caller file (TSV format, with 'sequence' column).")
    parser.add_argument("caller2", help="Path to the second variant caller file (TSV format, with 'EMPOP_Variant' column).")
    parser.add_argument("output_file", help="Path to save the merged output (TSV format).")

    args = parser.parse_args()

    df_merged = merge_variant_callers(args.caller1, args.caller2)
    df_merged.to_csv(args.output_file, sep="\t", index=False)
    # print(f"Merged variant table saved to: {args.output_file}")