import sys

def process_sam_line(line):
    # Keep header lines unchanged
    if line.startswith('@'):
        return line.strip()
    
    fields = line.strip().split('\t')
    
    # Ensure the line has at least the standard 11 fields for SAM format
    if len(fields) < 11:
        return None  # Skip malformed or empty lines

    cigar = fields[5]
    seq = fields[9]
    qual = fields[10]

    # Process only if soft-clipping exists
    if 'S' in cigar:
        import re
        # Remove soft-clipping from the start
        match_start = re.match(r'(\d+)S', cigar)
        if match_start:
            clip_len = int(match_start.group(1))
            seq = seq[clip_len:]
            qual = qual[clip_len:]
            cigar = cigar[len(match_start.group(0)):]  # Update CIGAR string

        # Remove soft-clipping from the end
        match_end = re.search(r'(\d+)S$', cigar)
        if match_end:
            clip_len = int(match_end.group(1))
            seq = seq[:-clip_len]
            qual = qual[:-clip_len]
            cigar = cigar[:-len(match_end.group(0))]  # Update CIGAR string

        # Update fields
        fields[5] = cigar
        fields[9] = seq
        fields[10] = qual

    # Skip reads where sequence and quality lengths don't match
    if len(seq) != len(qual):
        return None

    return '\t'.join(fields)

# Read from stdin and write to stdout
for line in sys.stdin:
    processed_line = process_sam_line(line)
    if processed_line:
        print(processed_line)

