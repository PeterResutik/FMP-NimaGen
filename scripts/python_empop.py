import pandas as pd
import argparse
from Bio import SeqIO

# IUPAC codes for heteroplasmy
IUPAC_CODES = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["A", "C"]): "M",
    frozenset(["G", "T"]): "K",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W"
}

def load_reference_genome(fasta_file):
    """Load reference genome from FASTA file."""
    record = SeqIO.read(fasta_file, "fasta")
    return str(record.seq)

def rightmost_position(reference, pos, variant, ref_base):
    """
    Finds the rightmost position where an insertion/deletion should be placed based on repeat structure.
    """
    pos = int(pos) - 1  # Convert to zero-based index
    ref_base = reference[pos] if pos < len(reference) else ref_base
    # print(reference[pos])
    
    # Move the insertion to the rightmost position in case of repeated bases
    while pos + 1 < len(reference) and reference[pos + 1] == variant:
        pos += 1

    return pos + 1  # Convert back to 1-based index

def rightmost_position_segment_ins(reference, pos, segment):
    """
    Finds the rightmost position where an insertion or deletion should be placed
    based on repeated sequences.
    """
    pos = int(pos) - 1  # Convert to zero-based index

    # Move position to the rightmost occurrence of the repeated segment
    while pos + len(segment) < len(reference) and reference[pos + 1: pos + 1 + len(segment)] == segment:
        pos += len(segment)
        print(pos)
        print(reference[pos + 1: pos + 1 + len(segment)])

    return pos+len(segment)  # Convert back to 1-based index

def rightmost_position_segment_del(reference, pos, segment):
    """
    Finds the rightmost position where an insertion or deletion should be placed
    based on repeated sequences.
    """
    pos = int(pos) - 1  # Convert to zero-based index

    # Move position to the rightmost occurrence of the repeated segment
    while pos + len(segment) < len(reference) and reference[pos + 1: pos + 1 + len(segment)] == segment:
        pos += len(segment)
        print(pos)
        print(reference[pos + 1: pos + 1 + len(segment)])

    return pos  # Convert back to 1-based index


def shift_insertion_right(reference, pos, inserted_segment):
    """
    Shifts the insertion right if the first base of inserted_segment matches 
    the next base after the rightmost position of the insertion.
    """
    rightmost_pos = rightmost_position_segment_ins(reference, pos, inserted_segment)

    for _ in range(len(inserted_segment) - 1):
        if rightmost_pos + len(inserted_segment) < len(reference) and reference[rightmost_pos + len(inserted_segment)] == inserted_segment[0]:
            # Shift right by moving first base to last position
            inserted_segment = inserted_segment[1:] + inserted_segment[0]
            rightmost_pos += 1
        else:
            break  # Stop shifting if no more matches

    return rightmost_pos, inserted_segment


def shift_deletion_right(reference, pos, deleted_segment):
    """
    Shifts the deletion right if the first base of deleted_segment matches 
    the next base after the rightmost position of the deletion.
    """
    rightmost_pos = rightmost_position_segment_del(reference, pos, deleted_segment)
    print(rightmost_pos)
    for _ in range(len(deleted_segment) - 1):
        if rightmost_pos + len(deleted_segment) < len(reference) and reference[rightmost_pos + len(deleted_segment)] == deleted_segment[0]:
            # Shift right by moving first base to last position
            deleted_segment = deleted_segment[1:] + deleted_segment[0]
            rightmost_pos += 1
        else:
            break  # Stop shifting if no more matches
    print(deleted_segment)
    print(rightmost_pos)
    return rightmost_pos, deleted_segment

def format_variant(row, reference):
    """
    Formats variants according to EMPOP standards.
    """
    pos, ref, var, var_level, var_type = row["Pos"], row["Ref"], row["Variant"], row["VariantLevel"], row["Type"]
    pos = int(pos)  # Convert position to integer

    if var_type == 2:  # SNPs
        if var_level >= 0.925:
            return f"{ref}{pos}{var}"
        else:
            iupac_code = IUPAC_CODES.get(frozenset([ref, var]), f"{ref}/{var}")
            return f"{ref}{pos}{iupac_code}"

    elif var_type == 3:  # Indels
        if len(ref) + 1 == len(var):  # Single-base Insertion
            
            inserted_base = var[len(ref):]  # Extract inserted base
            rightmost_pos = rightmost_position(reference, pos, inserted_base, ref[-1])  # Find rightmost position
            if var_level >= 0.925:
                return f"-{rightmost_pos}.1{inserted_base}"  # Example: T315.1C
            else: 
                return f"-{rightmost_pos}.1{inserted_base.lower()}" 

        elif len(ref) - 1 == len(var):  # Single-base Deletion
            deleted_base = ref[len(var):]  # Extract deleted base
            rightmost_pos = rightmost_position(reference, pos, deleted_base, ref[-1])  # Find rightmost position
            return f"{deleted_base}{rightmost_pos}-"  # Example: C524-

        elif len(ref) < len(var):  # Multi-base Insertion
            inserted_segment = var[len(ref):]  # Extract inserted bases
            rightmost_pos, adjusted_segment = shift_insertion_right(reference, pos, inserted_segment)  # Find rightmost repeat
            if var_level >= 0.925:
                return " ".join([f"-{rightmost_pos}.{i+1}{adjusted_segment[i]}" for i in range(len(inserted_segment))])
            else: 
                return " ".join([f"-{rightmost_pos}.{i+1}{adjusted_segment[i].lower()}" for i in range(len(inserted_segment))])
        
        elif len(ref) > len(var):  # Multi-base Deletion
            deleted_segment = ref[len(var):]  # Extract deleted bases
            print(deleted_segment)
            # rightmost_pos = rightmost_position_segment_del(reference, pos, deleted_segment)  # Find rightmost repeat
            rightmost_pos, adjusted_segment = shift_deletion_right(reference, pos, deleted_segment)  # Find rightmost repeat
            return " ".join([f"{adjusted_segment[i]}{rightmost_pos + i}-" for i in range(len(adjusted_segment))])
    
    # var_type == 3:  # Indels
    #     if len(ref) < len(var):  # Insertion
    #         insertion_base = var[len(ref):]
    #         print(insertion_base)
            
    #         rightmost_pos = rightmost_position(reference, pos, insertion_base, ref[-1])
    #         # print(rightmost_pos)
    #         return " ".join([f"{ref}{rightmost_pos}.{i+1}{insertion_base[i]}" for i in range(len(insertion_base))])

    #     elif len(ref) > len(var):  # Deletion
    #         deleted_bases = ref[len(var):]
    #         rightmost_pos = rightmost_position(reference, pos, deleted_bases, ref[-1])
    #         return " ".join([f"{deleted_bases[i]}{rightmost_pos + i}-" for i in range(len(deleted_bases))])

    return "N/A"  # If type does not match expected cases
    
def process_variants(input_file, output_file, ref_fasta):
    """
    Reads variant file, processes variants into EMPOP format, and saves to output file.
    """
    # Load reference genome
    reference = load_reference_genome(ref_fasta)

    # Load variant data
    df = pd.read_csv(input_file, sep="\t")

    # Ensure VariantLevel is float
    df["VariantLevel"] = pd.to_numeric(df["VariantLevel"], errors="coerce").fillna(1.0)

    # Process variants
    df["EMPOP_Variant"] = df.apply(lambda row: format_variant(row, reference), axis=1)

    # Save updated file
    df.to_csv(output_file, sep="\t", index=False)

    print(f"Processed file saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Format mitochondrial variants for EMPOP.")
    parser.add_argument("input_file", help="Path to the input variant file (TSV format).")
    parser.add_argument("output_file", help="Path to save the formatted output file (TSV format).")
    parser.add_argument("ref_fasta", help="Path to the mitochondrial reference genome (FASTA format).")

    args = parser.parse_args()

    process_variants(args.input_file, args.output_file, args.ref_fasta)
