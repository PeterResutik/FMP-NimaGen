#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import os

def main():
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Plot mtDNA Amplicon Read Depth Comparison")
    parser.add_argument("input_file1", help="Path to the first read depth file (TSV format)")
    parser.add_argument("input_file2", help="Path to the second read depth file (TSV format)")
    parser.add_argument("output_file", help="Path to save the output plot (e.g., read_depth_comparison.png)")
    args = parser.parse_args()

    try:
        # Check file existence
        if not os.path.isfile(args.input_file1) or not os.path.isfile(args.input_file2):
            raise FileNotFoundError("One or both input files not found.")

        # Load the read depth data
        columns = ["Chromosome", "Position", "ReadDepth"]
        data1 = pd.read_csv(args.input_file1, sep="\t", header=None, names=columns)
        data2 = pd.read_csv(args.input_file2, sep="\t", header=None, names=columns)

        # Merge both dataframes on Position
        merged_data = pd.merge(data1, data2, on="Position", suffixes=("_file1", "_file2"))

        # Sort by position
        merged_data = merged_data.sort_values(by="Position")

        # Define read depth drop (safe division)
        merged_data["ReadDepth_Drop"] = merged_data.apply(
            lambda row: ((row["ReadDepth_file1"] - row["ReadDepth_file2"]) / row["ReadDepth_file1"] * 100)
            if row["ReadDepth_file1"] > 0 else 0,
            axis=1
        )

        # Set color mapping (not used in plot, but retained)
        colors = [
            "red" if drop > 50 else "orange" if drop > 25 else "yellow" if drop > 10 else "blue"
            for drop in merged_data["ReadDepth_Drop"]
        ]

        # Plot as a stacked bar plot
        plt.figure(figsize=(12, 6))
        plt.bar(merged_data["Position"], merged_data["ReadDepth_file1"], color="blue", label="Read Depth w/ NUMTs", width=50)
        plt.bar(merged_data["Position"], merged_data["ReadDepth_file2"], color="red", label="Read Depth w/o NUMTs", width=50, alpha=0.7)

        # Log scale for read depth
        plt.yscale("log")

        # Labels and title
        plt.xlabel("Amplicon Middle Position (chrM)")
        plt.ylabel("Read Depth (log scale)")
        plt.title("Comparison of mtDNA Amplicon Read Depth")

        # Add horizontal dashed lines at read depth 10 and 20
        plt.axhline(y=10, color="red", linestyle="--", linewidth=1, label="Read Depth = 10")
        plt.axhline(y=20, color="gray", linestyle="--", linewidth=1, label="Read Depth = 20")

        # Grid settings
        plt.grid(axis="y", linestyle="--", linewidth=0.5)

        # Add x-axis ticks at every amplicon position
        positions = merged_data["Position"].tolist()
        plt.xticks(positions, labels=[str(pos) for pos in positions], rotation=90, fontsize=5)

        # Legend and save
        plt.legend()
        plt.savefig(args.output_file, dpi=300, bbox_inches="tight")

        print(f"Plot saved as {args.output_file}")

    except FileNotFoundError as e:
        print(f"File error: {e}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("One of the input files is empty or improperly formatted.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
