import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Command-line argument parsing
parser = argparse.ArgumentParser(description="Plot mtDNA Amplicon Coverage Comparison")
parser.add_argument("input_file1", help="Path to the first coverage file (TSV format)")
parser.add_argument("input_file2", help="Path to the second coverage file (TSV format)")
parser.add_argument("output_file", help="Path to save the output plot (e.g., coverage_comparison.png)")
args = parser.parse_args()

# Load the coverage data
columns = ["Chromosome", "Position", "Coverage"]
data1 = pd.read_csv(args.input_file1, sep="\t", header=None, names=columns)
data2 = pd.read_csv(args.input_file2, sep="\t", header=None, names=columns)

# Merge both dataframes on Position
merged_data = pd.merge(data1, data2, on="Position", suffixes=("_file1", "_file2"))

# Sort by position
merged_data = merged_data.sort_values(by="Position")

# Define colors based on the percentage drop in coverage
merged_data["Coverage_Drop"] = (merged_data["Coverage_file1"] - merged_data["Coverage_file2"]) / merged_data["Coverage_file1"] * 100

# Set color mapping
colors = [
    "red" if drop > 50 else "orange" if drop > 25 else "yellow" if drop > 10 else "blue"
    for drop in merged_data["Coverage_Drop"]
]

# Plot as a stacked bar plot
plt.figure(figsize=(12, 6))
plt.bar(merged_data["Position"], merged_data["Coverage_file1"], color="blue", label="Coverage in File 1", width=50)
plt.bar(merged_data["Position"], merged_data["Coverage_file2"], color="red", label="Coverage in File 2", width=50, alpha=0.7)

# Log scale for coverage
plt.yscale("log")

# Labels and title
plt.xlabel("Amplicon Middle Position (chrM)")
plt.ylabel("Coverage (log scale)")
plt.title("Comparison of mtDNA Amplicon Coverage")

# Add horizontal dashed lines at coverage 10 and 20
plt.axhline(y=10, color="red", linestyle="--", linewidth=1, label="Coverage = 10")
plt.axhline(y=20, color="gray", linestyle="--", linewidth=1, label="Coverage = 20")

# Grid settings
plt.grid(axis="y", linestyle="--", linewidth=0.5)

# Add x-axis ticks at every amplicon position
positions = merged_data["Position"].tolist()
plt.xticks(positions, labels=[str(pos) for pos in positions], rotation=90, fontsize=5)

# Legend and save
plt.legend()
plt.savefig(args.output_file, dpi=300, bbox_inches="tight")

print(f"Plot saved as {args.output_file}")