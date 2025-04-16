#!/usr/bin/env python3

import sys
import re

# Precompiled regex patterns
SOFT_CLIP_START = re.compile(r'^(\d+)S')
SOFT_CLIP_END = re.compile(r'(\d+)S$')

def process_sam_line(line):
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

        return '\t'.join(fields)
    
    except Exception as e:
        sys.stderr.write(f"Error processing line:\n{line}\n{e}\n")
        return None

def main():
    if len(sys.argv) != 3:
        print("Usage: python remove_scb.py input.sam output.sam", file=sys.stderr)
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line_number, line in enumerate(infile, start=1):
                processed = process_sam_line(line)
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
