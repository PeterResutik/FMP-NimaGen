#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import re
import argparse
import sys

import frames
import notation

IUPAC_CODES = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["A", "C"]): "M",
    frozenset(["G", "T"]): "K",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W"
}
IUPAC_BASES = {code: set(bases) for bases, code in IUPAC_CODES.items()}

def claim_key(label):
    """Same claim regardless of major/minor spelling (T16189C/T16189Y,
    -309.1C/-309.1c, A523-/A523a); different alleles at one position
    (C756M vs C756Y) get different keys."""
    label = str(label)
    m = re.match(r"^([ACGT])(\d+)([A-Za-z-])$", label)
    if m:
        ref, pos, code = m.groups()
        if code == "-" or code == ref.lower():
            return f"{ref}{pos}-"
        others = IUPAC_BASES.get(code, set()) - {ref}
        if len(others) == 1:
            return f"{ref}{pos}{others.pop()}"
    return label.upper()

# rCRS_NimaGen.fasta = 16569bp rCRS + an appended copy of chrM:1-53, so the
# origin-spanning amplicon maps linearly. A POS past 16569 is the same base
# as POS-16569.
TRUE_MT_LENGTH = 16569

def wrap_circular_position(pos, true_length=TRUE_MT_LENGTH):
    return pos - true_length if pos > true_length else pos

def pool_origin_overlap(df):
    """Merge a variant seen through both the appended copy and its true
    coordinate (two amplicons, independent reads) into one record, as a
    single pileup would: AD summed, AF recalculated as alt / total reads."""
    if df.empty:
        return df
    df = df.copy()
    df["_wrapped_pos"] = df["Pos"].astype(int).apply(wrap_circular_position)
    df["_depth"] = df["Coverage"].astype(str).apply(lambda ad: sum(float(x) for x in ad.split(",")))
    pooled = []
    for _, group in df.groupby(["_wrapped_pos", "Ref", "Variant"], sort=False):
        if len(group) == 1:
            pooled.append(group.iloc[0])
            continue
        row = group.loc[group["_depth"].idxmax()].copy()
        ad_sums = [sum(x) for x in zip(*group["Coverage"].astype(str).apply(lambda ad: [float(v) for v in ad.split(",")]))]
        total_depth = sum(ad_sums)
        row["Coverage"] = ",".join(f"{v:.4g}" for v in ad_sums)
        if total_depth:
            row["VariantLevel"] = round(ad_sums[1] / total_depth, 4)
        row["Pos"] = row["_wrapped_pos"]
        pooled.append(row)
    return pd.DataFrame(pooled).drop(columns=["_wrapped_pos", "_depth"]).infer_objects().reset_index(drop=True)

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

def lh_bounds(threshold):
    """Symmetric floor/ceiling around a single length-heteroplasmy threshold,
    e.g. threshold=0.10 -> (0.10, 0.90). Accepts either side (0.10 or 0.90)
    and always returns (floor, ceiling) with floor <= ceiling."""
    return min(threshold, 1 - threshold), max(threshold, 1 - threshold)

def shift_deletion_right(reference, pos, segment):
    pos = rightmost_repeat_position(reference, pos, segment) - len(segment) + 1
    for _ in range(len(segment) - 1):
        next_seq = reference[pos + len(segment) - 1: pos + len(segment)]
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
        label_pos = wrap_circular_position(sub_pos)
        if var_level >= 1 - min_variant_frequency:
            formatted.append(f"{r}{label_pos}{v}")
        else:
            code = IUPAC_CODES.get(frozenset([r, v]), f"{r}/{v}")
            formatted.append(f"{r}{label_pos}{code}")

    # Return after the loop finishes
    if var_level >= 1 - min_variant_frequency:
        return " ".join(formatted), "SNP"
    else:
        return " ".join(formatted), "PHP"

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
    deleted_segment = ref[len(var):]
    if pos == 16188 and reference[pos]=="C":
        deleted_segment = "".join(reference[pos:pos+len(var)])
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
    df["_claim"] = df["MUTECT2"].apply(claim_key)

    # One row per claim: identical claims (e.g. two alleles of one insertion landing
    # on the same shifted label) are summed; different alleles at one position are not
    group_keys = ["_claim"]
    numeric_agg = {
        "VariantLevel": "sum",
        "Coverage": add_comma_separated_numbers,
        "MeanBaseQuality": "first"
    }
    other_cols = [col for col in df.columns if col not in numeric_agg and col not in group_keys]
    full_agg = {**numeric_agg, **{col: "first" for col in other_cols}}

    grouped = df.groupby(group_keys, as_index=False).agg(full_agg)

    # Optional: sort and drop temp column
    grouped = grouped.sort_values("variant_float_pos").drop(columns=["variant_float_pos", "_claim"])


    # def extract_position(variant):
    #     match = re.search(r"(\d+\.?\d*)", str(variant))
    #     return float(match.group(1)) if match else float('inf')

    def extract_position(variant):
        match = re.search(r"(\d+\.?\d*)", str(variant))
        if match:
            return int(float(match.group(1)))
        return float('inf')

    grouped["position"] = grouped["MUTECT2"].apply(extract_position)
    grouped = grouped.sort_values(by="position").drop(columns=["position"])

    _, lh_ceiling = lh_bounds(length_heteroplasmy_threshold)

    def correct_length_het_case(row):
        if "." in row["MUTECT2"] and row["VariantLevel"] >= lh_ceiling:
            return row["MUTECT2"][:-1] + row["MUTECT2"][-1].upper()
        return row["MUTECT2"]

    grouped["MUTECT2"] = grouped.apply(correct_length_het_case, axis=1)
    return grouped

def label_type(label):
    return "INS" if label.startswith("-") else "DEL" if label.endswith("-") else "SNP"


def replace_labels(df, old, new):
    """Rows old {index: label} rewritten as the labels new: a row whose label is among the
    new ones keeps its frequency (and takes that label), the other rows go, and each new
    label left over gets the lowest frequency of the rows that went. Type follows the label."""
    staying = {i: label for i, label in old.items() if label in new}
    going = [i for i in old if i not in staying]
    added = [label for label in new if label not in staying.values()]
    df = df.copy()
    for i, label in staying.items():
        df.loc[i, "MUTECT2"], df.loc[i, "Type"] = label, label_type(label)
    rows = []
    if added:
        lowest = df.loc[going or list(old)].sort_values("VariantLevel").iloc[0].to_dict()
        rows = [{**lowest, "MUTECT2": label, "Type": label_type(label)} for label in added]
    return pd.concat([df.drop(index=going), pd.DataFrame(rows)], ignore_index=True)


def respell_rows(df, reference):
    """Groups of major labels outside the frame regions written again by the general
    rule (notation.respell_majors), so both callers spell a molecule alike."""
    protected = [(r.first, r.last) for r in frames.REGIONS.values()]
    for old, new in notation.respell_majors(df["MUTECT2"], reference, protected):
        df = replace_labels(df, {i: label for i, label in df["MUTECT2"].items() if label in old}, new)
    return df


def frame_rows(df, reference, frame):
    """Mutect2's major calls in each frame region and within frames.FLANK bases of it,
    written as FDSTOOLS' molecules are: rebuilt into one molecule and placed by the same
    alignment, the region in the chosen frame and its edges by the general rule. Minor
    calls stay as Mutect2 wrote them."""
    for region in frames.REGIONS.values():
        lo, hi = region.first - frames.FLANK, region.last + frames.FLANK
        majors = {}
        for idx, row in df.iterrows():
            label = str(row["MUTECT2"])
            if not lo <= notation.label_position(label)[0] <= hi:
                continue
            if notation.is_major(label):
                majors[idx] = label
        if not majors:
            continue
        molecule = notation.apply_labels(majors.values(), reference, lo, hi)
        bases, flank = frames.region_and_flank(molecule, lo, hi, region, reference)
        new = frames.labels(region, bases, frame, reference) + flank
        if sorted(new) != sorted(df.loc[list(majors), "MUTECT2"]):
            df = replace_labels(df, majors, new)
    return df


def main():
    parser = argparse.ArgumentParser(description="Process mitochondrial variants into EMPOP format.")
    parser.add_argument("input_file", help="Input TSV file with variants")
    parser.add_argument("output_file", help="Output TSV file")
    parser.add_argument("reference_fasta", help="Reference genome in FASTA format")
    parser.add_argument("--min_vf", type=float, default=5.0, help="Minor allele frequency threshold")
    parser.add_argument("--lh_thresh", type=float, default=10.0, help="Symmetric length-heteroplasmy threshold (floor and, via 1-threshold, ceiling), e.g. 10.0 -> report only 10-90%%, lowercase in between, major above 90%%")
    parser.add_argument("--frame", choices=["separate", "shared"], default="separate", help="Frame for 16180-16193 and 300-315, as for FDSTOOLS")

    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_file, sep="\t")
        df = pool_origin_overlap(df)
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

        # Drop length variants below the symmetric length-heteroplasmy floor entirely
        df = df[df["Type"] != "BELOW_LH_FLOOR"].reset_index(drop=True)

        df = finalize_output_table(df, args.lh_thresh/100)
        df = respell_rows(df, load_reference(args.reference_fasta))
        df = frame_rows(df, load_reference(args.reference_fasta), args.frame)

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
