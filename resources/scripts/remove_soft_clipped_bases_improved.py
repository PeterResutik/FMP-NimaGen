#!/usr/bin/env python3

import sys
import re

# Precompiled regex patterns
SOFT_CLIP_START = re.compile(r'^(\d+)S')
SOFT_CLIP_END = re.compile(r'(\d+)S$')
CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')

# Bases a read reaches past a primer boundary: read-through into the kit's tail that happens to
# match the reference, so bwa aligns it and soft-clip removal keeps it. It is a few bases at most.
READ_THROUGH = 5
IUPAC = {"R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]", "S": "[CG]", "W": "[AT]",
         "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}


def fasta_sequences(path):
    seqs, current = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
            current = []
        elif line:
            current.append(line.upper().strip("^$"))
    if current:
        seqs.append("".join(current))
    return seqs


def primer_boundaries(reference, left_primers, right_primers):
    """Reference positions where a product can begin (left-primer starts) and end (right-primer ends);
    the right primers are given as on the forward strand."""
    ref = "".join(fasta_sequences(reference))
    def matches(primer):
        return re.finditer("(?=" + "".join(IUPAC.get(b, b) for b in primer) + ")", ref)
    lefts = sorted({m.start() + 1 for p in fasta_sequences(left_primers) for m in matches(p)})
    rights = sorted({m.start() + len(p) for p in fasta_sequences(right_primers) for m in matches(p)})
    return lefts, rights


def trim_cigar(ops, ref_bases):
    """Remove ref_bases reference bases from the start of ops (a list of (length, op));
    returns the remaining ops and the query and reference bases removed."""
    ops, q_cut, r_cut = list(ops), 0, 0
    while ops and (r_cut < ref_bases or ops[0][1] in "ID"):
        n, op = ops.pop(0)
        if op in "M=X":
            k = min(n, ref_bases - r_cut)
            if k == 0:
                ops.insert(0, (n, op))
                break
            q_cut, r_cut = q_cut + k, r_cut + k
            if n > k:
                ops.insert(0, (n - k, op))
        elif op == "I":
            q_cut += n
        elif op in "DN":
            r_cut += n
    return ops, q_cut, r_cut


def trim_read_through(fields, lefts, rights, window=READ_THROUGH):
    """Cut a mapped read back to the nearest primer boundary it reaches past by at most window bases.
    A product spanning two amplicons begins and ends at primer boundaries too, so it stays whole."""
    if int(fields[1]) & (4 | 256 | 2048) or fields[5] in ("*", ""):
        return fields
    ops = [(int(n), op) for n, op in CIGAR.findall(fields[5])]
    if not ops or ops[0][1] in "SH" or ops[-1][1] in "SH":
        return fields
    start = int(fields[3])
    end = start + sum(n for n, op in ops if op in "MDN=X") - 1
    seq, qual = fields[9], fields[10]
    left = next((b for b in lefts if 0 < b - start <= window), None)
    if left is not None:
        ops, q_cut, r_cut = trim_cigar(ops, left - start)
        seq, qual, start = seq[q_cut:], qual[q_cut:], start + r_cut
    right = next((b for b in rights if 0 < end - b <= window), None)
    if right is not None:
        rev, q_cut, _ = trim_cigar(ops[::-1], end - right)
        ops = rev[::-1]
        seq, qual = seq[:len(seq) - q_cut], qual[:len(qual) - q_cut]
    if not ops or not seq:
        return fields
    fields[3], fields[5], fields[9], fields[10] = str(start), "".join(f"{n}{op}" for n, op in ops), seq, qual
    return fields

def process_sam_line(line, boundaries=None):
    try:
        if line.startswith('@'):
            return line.rstrip('\n')

        fields = line.rstrip('\n').split('\t')
        if len(fields) < 11:
            return None

        cigar = fields[5]
        seq = fields[9]
        qual = fields[10]

        if 'S' in cigar:
            match_start = SOFT_CLIP_START.match(cigar)
            if match_start:
                clip_len = int(match_start.group(1))
                seq = seq[clip_len:]
                qual = qual[clip_len:]
                cigar = cigar[len(match_start.group(0)):]

            match_end = SOFT_CLIP_END.search(cigar)
            if match_end:
                clip_len = int(match_end.group(1))
                seq = seq[:-clip_len] if clip_len < len(seq) else ''
                qual = qual[:-clip_len] if clip_len < len(qual) else ''
                cigar = cigar[:-len(match_end.group(0))]

            fields[5] = cigar
            fields[9] = seq
            fields[10] = qual

        if len(seq) != len(qual):
            return None

        if boundaries:
            fields = trim_read_through(fields, *boundaries)

        return '\t'.join(fields)
    
    except Exception as e:
        sys.stderr.write(f"Error processing line:\n{line}\n{e}\n")
        return None

def main():
    if len(sys.argv) not in (3, 6):
        print("Usage: python remove_scb.py input.sam output.sam [reference.fasta left_primers.fasta right_primers_rc.fasta]",
              file=sys.stderr)
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]
    boundaries = primer_boundaries(*sys.argv[3:6]) if len(sys.argv) == 6 else None

    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line_number, line in enumerate(infile, start=1):
                processed = process_sam_line(line, boundaries)
                if processed:
                    outfile.write(processed + '\n')
    except FileNotFoundError:
        print(f"Error: File not found '{input_path}'", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"I/O error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
