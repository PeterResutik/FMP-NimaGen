#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import re
import argparse
import sys

IUPAC_CODES = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["A", "C"]): "M",
    frozenset(["G", "T"]): "K",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W"
}

def load_reference(fasta_path):
    record = SeqIO.read(fasta_path, "fasta")
    return list(str(record.seq))

def rightmost_repeat_position(reference, pos, segment):
    pos -= 1
    while (
        pos + len(segment) < len(reference) and
        "".join(reference[pos + 1: pos + 1 + len(segment)]) == segment
    ):
        pos += len(segment)
    return pos + 1

def shift_insertion_right(reference, pos, segment):
    pos = rightmost_repeat_position(reference, pos, segment)
    for _ in range(len(segment) - 1):
        next_seq = reference[pos: pos + 1]
        if "".join(next_seq) == segment[0]:
            segment = segment[1:] + segment[0]
            pos += 1
        else:
            break
    return pos, segment

def shift_deletion_right(reference, pos, segment):
    pos = rightmost_repeat_position(reference, pos, segment) - len(segment) + 1
    next_seq = reference[pos+1: pos + 2]
    for _ in range(len(segment) - 1):
        if "".join(next_seq) == segment[0]:
            segment = segment[1:] + segment[0]
            pos += 1
        else:
            break
    return pos, segment

def apply_snp(pos, ref, var, var_level, reference, min_variant_frequency):
    formatted = []
    for i, (r, v) in enumerate(zip(ref, var)):
        sub_pos = pos + i
        reference[sub_pos - 1] = v
        if var_level >= 1 - min_variant_frequency:
            formatted.append(f"{r}{sub_pos}{v}")
        else:
            code = IUPAC_CODES.get(frozenset([r, v]), f"{r}/{v}")
            formatted.append(f"{r}{sub_pos}{code}")

    # Return after the loop finishes
    if var_level >= 1 - min_variant_frequency:
        return " ".join(formatted), "SNP"
    else:
        return " ".join(formatted), "PHP"

def apply_insertion(pos, ref, var, var_level, reference, length_heteroplasmy_threshold):
    inserted_segment = var[len(ref):]
    pos, segment = shift_insertion_right(reference, pos, inserted_segment)
    variant_parts = [
        f"-{pos}.{i+1}{(b if var_level >= length_heteroplasmy_threshold else b.lower())}"
        for i, b in enumerate(segment)
    ]
    updated_type = "INS" if var_level >= length_heteroplasmy_threshold else "LHP"
    return " ".join(variant_parts), updated_type

def apply_deletion(pos, ref, var, var_level, reference, length_heteroplasmy_threshold):
    deleted_segment = ref[len(var):]
    if pos == 16188 and reference[pos]=="C":
        deleted_segment = "".join(reference[pos:pos+len(var)])
    pos, segment = shift_deletion_right(reference, pos, deleted_segment)

    is_major = var_level >= length_heteroplasmy_threshold
    variant_parts = []
    for i, base in enumerate(segment):
        position = pos + i
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

    # group_keys = ["MUTECT2"]
    # numeric_agg = {
    #     "VariantLevel": "sum",
    #     "Coverage": add_comma_separated_numbers,
    #     "MeanBaseQuality": "first"
    # }
    # other_cols = [col for col in df.columns if col not in numeric_agg and col not in group_keys]
    # full_agg = {**numeric_agg, **{col: "first" for col in other_cols}}

    # grouped = df.groupby(group_keys, as_index=False).agg(full_agg)



    # Add column for numeric variant position
    df["variant_float_pos"] = df["MUTECT2"].apply(extract_float_position)

    group_keys = ["variant_float_pos"]  # <-- changed here
    numeric_agg = {
        "VariantLevel": "sum",
        "Coverage": add_comma_separated_numbers,
        "MeanBaseQuality": "first"
    }
    other_cols = [col for col in df.columns if col not in numeric_agg and col not in group_keys]
    full_agg = {**numeric_agg, **{col: "first" for col in other_cols}}

    grouped = df.groupby(group_keys, as_index=False).agg(full_agg)

    # Optional: sort and drop temp column
    grouped = grouped.sort_values("variant_float_pos").drop(columns=["variant_float_pos"])


    # def extract_position(variant):
    #     match = re.search(r"(\d+\.?\d*)", str(variant))
    #     return float(match.group(1)) if match else float('inf')

    def extract_position(variant):
        match = re.search(r"(\d+\.?\d*)", str(variant))
        if match:
            pos = int(float(match.group(1)))
            if 16570 <= pos <= 16587:
                return pos - 16569  # map 16570–16587 → 1–18
            return pos
        return float('inf')

    grouped["position"] = grouped["MUTECT2"].apply(extract_position)
    grouped = grouped.sort_values(by="position").drop(columns=["position"])

    def correct_length_het_case(row):
        if "." in row["MUTECT2"] and row["VariantLevel"] >= length_heteroplasmy_threshold:
            return row["MUTECT2"][:-1] + row["MUTECT2"][-1].upper()
        return row["MUTECT2"]

    grouped["MUTECT2"] = grouped.apply(correct_length_het_case, axis=1)
    return grouped

def main():
    parser = argparse.ArgumentParser(description="Process mitochondrial variants into EMPOP format.")
    parser.add_argument("input_file", help="Input TSV file with variants")
    parser.add_argument("output_file", help="Output TSV file")
    parser.add_argument("reference_fasta", help="Reference genome in FASTA format")
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minor allele frequency threshold")
    parser.add_argument("--lh_thresh", type=float, default=90.0, help="Length heteroplasmy frequency threshold")

    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_file, sep="\t")
        reference = load_reference(args.reference_fasta)

        results = []
        types = []
        for idx in reversed(df.index):
            row = df.loc[idx]
            variant, updated_type = format_variant(row, reference, args.min_vf/100, args.lh_thresh/100)
            results.append(variant)
            types.append(updated_type)

        df["MUTECT2"] = results[::-1]
        df["Type"] = types[::-1]

        df = finalize_output_table(df, args.lh_thresh/100)

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
