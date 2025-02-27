import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Command-line argument parsing
parser = argparse.ArgumentParser(description="Plot mtDNA Amplicon Coverage on a Log Scale")
parser.add_argument("input_file", help="Path to the coverage file (TSV format)")
parser.add_argument("output_file", help="Path to save the output plot (e.g., coverage_plot.png)")
args = parser.parse_args()

# Load the coverage data
data = pd.read_csv(args.input_file, sep="\t", header=None, names=["Chromosome", "Position", "Coverage"])

# Sort by position for better visualization
data = data.sort_values(by="Position")

# Ensure all amplicon positions have x-axis ticks, even if missing in the dataset
positions = data["Position"].tolist()

# Define threshold for low coverage
low_coverage_threshold = 100
colors = ["red" if cov < low_coverage_threshold else "blue" for cov in data["Coverage"]]

# Plot as a bar plot
plt.figure(figsize=(12, 6))
plt.bar(data["Position"], data["Coverage"], color=colors, width=50, label="Coverage")

# Logarithmic scale for coverage
plt.yscale("log")

# Add horizontal dashed lines at coverage 10 and 20
plt.axhline(y=10, color="black", linestyle="--", linewidth=1, label="Coverage = 10")
plt.axhline(y=20, color="gray", linestyle="--", linewidth=1, label="Coverage = 20")

# Labels and title
plt.xlabel("Amplicon Middle Position (chrM)")
plt.ylabel("Coverage (log scale)")
plt.title("mtDNA Amplicon Coverage")

# Remove vertical grid lines while keeping horizontal ones
plt.grid(axis="y", linestyle="--", linewidth=0.5)

# Add x-axis ticks at every amplicon position (so 0-coverage bars are recognizable)
plt.xticks(positions, labels=[str(pos) for pos in positions], rotation=90, fontsize=5)


# Adjust x-ticks
plt.xticks(rotation=45, fontsize=10)
plt.legend()

# Save the plot to a file
plt.savefig(args.output_file, dpi=300, bbox_inches="tight")

# print(f"Plot saved as {args.output_file}")
