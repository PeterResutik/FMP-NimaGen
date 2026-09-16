#!/usr/bin/env python3
"""Direct read-level evidence from a BAM, for loci where variant-caller
allele bookkeeping is unreliable.

Motivation: in homopolymer/repeat-adjacent regions, both Mutect2's
VCF allele frequencies and FDSTOOLS' haplotype counts can disagree
sharply with what the reads actually show - fragmented multiallelic
splits whose AFs don't sum sensibly, local-reassembly artifacts that
invent indel lengths present in no read, and (on the FDSTOOLS side)
del+sub merges that fold distinct events together. This module goes
back to the reads and counts, with no quality or mapping-quality
filtering, since at these loci the reference-matching reads are exactly
the ones that get preferentially filtered out (samtools mpileup's
default -Q 13 silently removed most of them in one case, inverting the
apparent picture).

Everything here is a plain CIGAR walk - no pysam dependency, no
reliance on samtools' pileup filtering semantics.
"""

import argparse
import re
import subprocess
import sys

CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def _iter_reads(bam_path, region):
    """Yield (pos, cigar, seq) for each alignment overlapping region.

    Uses `samtools view` with no filtering flags, so nothing is excluded
    on quality grounds - that filtering is precisely what misleads at
    these loci.
    """
    proc = subprocess.run(
        ["samtools", "view", bam_path, region],
        capture_output=True, text=True, check=True,
    )
    for line in proc.stdout.splitlines():
        if not line or line.startswith("@"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        flag = int(fields[1])
        if flag & 0x100 or flag & 0x800:  # secondary / supplementary
            continue
        yield int(fields[3]), fields[5], fields[9]


def _walk(pos, cigar):
    """Yield (ref_pos, read_idx, op) for each consumed position.

    ref_pos is None for inserted bases (they consume no reference);
    read_idx is None for deleted positions (they consume no read base).
    """
    ref_pos = pos
    read_idx = 0
    for length, op in CIGAR_RE.findall(cigar):
        length = int(length)
        if op in ("M", "=", "X"):
            for i in range(length):
                yield ref_pos + i, read_idx + i, op
            ref_pos += length
            read_idx += length
        elif op == "I":
            for i in range(length):
                yield None, read_idx + i, op
            read_idx += length
        elif op in ("D", "N"):
            for i in range(length):
                yield ref_pos + i, None, op
            ref_pos += length
        elif op == "S":
            read_idx += length
        elif op == "H":
            pass


def base_composition(bam_path, chrom, position):
    """Observed base counts at a 1-based reference position.

    Returns {"A": n, "C": n, "G": n, "T": n, "deleted": n, "total": n},
    where "deleted" counts reads spanning the position with a deletion
    (they genuinely have no base there, which is different from showing
    the reference base).
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "deleted": 0, "total": 0}
    region = f"{chrom}:{position}-{position}"
    for pos, cigar, seq in _iter_reads(bam_path, region):
        for ref_pos, read_idx, op in _walk(pos, cigar):
            if ref_pos != position:
                continue
            counts["total"] += 1
            if op in ("D", "N"):
                counts["deleted"] += 1
            elif read_idx is not None:
                base = seq[read_idx].upper()
                if base in counts:
                    counts[base] += 1
            break
    return counts


def net_length_change(bam_path, chrom, start, end):
    """Distribution of net inserted/deleted bases per read within a window.

    Deliberately position-agnostic inside the window: in a homopolymer
    the exact CIGAR offset of an indel is arbitrary, but the net length
    change is not. Returns {net_change: read_count}, e.g. {0: 34, -1: 17,
    1: 13, -2: 7, 2: 2} meaning 34 reads span the window at reference
    length, 17 are one base short, 13 one base long, and so on.

    Only reads fully spanning the window are counted, so partial
    coverage can't masquerade as a length change.
    """
    dist = {}
    region = f"{chrom}:{start}-{end}"
    for pos, cigar, seq in _iter_reads(bam_path, region):
        covered_start = False
        covered_end = False
        net = 0
        for ref_pos, read_idx, op in _walk(pos, cigar):
            if op == "I":
                # An insertion sits between reference positions; count it
                # if the surrounding alignment is inside the window.
                if start <= _last_ref_pos(pos, cigar, read_idx) <= end:
                    net += 1
                continue
            if ref_pos is None:
                continue
            if ref_pos == start:
                covered_start = True
            if ref_pos == end:
                covered_end = True
            if start <= ref_pos <= end and op in ("D", "N"):
                net -= 1
        if covered_start and covered_end:
            dist[net] = dist.get(net, 0) + 1
    return dict(sorted(dist.items()))


def _last_ref_pos(pos, cigar, read_idx_at_insertion):
    """Reference position immediately preceding the read base at read_idx."""
    ref_pos = pos
    idx = 0
    last = pos
    for length, op in CIGAR_RE.findall(cigar):
        length = int(length)
        if op in ("M", "=", "X"):
            if idx <= read_idx_at_insertion < idx + length:
                return ref_pos + (read_idx_at_insertion - idx)
            ref_pos += length
            idx += length
            last = ref_pos - 1
        elif op == "I":
            if idx <= read_idx_at_insertion < idx + length:
                return last
            idx += length
        elif op in ("D", "N"):
            ref_pos += length
            last = ref_pos - 1
        elif op == "S":
            idx += length
    return last


def main():
    parser = argparse.ArgumentParser(
        description="Direct read-level evidence at a locus (no quality filtering)."
    )
    parser.add_argument("bam", help="Indexed BAM file")
    parser.add_argument("--chrom", default="chrM")
    parser.add_argument("--position", type=int, help="Report base composition at this 1-based position")
    parser.add_argument("--window", help="Report net length change over START-END, e.g. 16180-16193")
    args = parser.parse_args()

    if args.position:
        counts = base_composition(args.bam, args.chrom, args.position)
        called = {b: counts[b] for b in "ACGT" if counts[b]}
        total_called = sum(called.values())
        print(f"{args.chrom}:{args.position}  total_spanning={counts['total']}  deleted={counts['deleted']}")
        for base, n in sorted(called.items(), key=lambda kv: -kv[1]):
            pct = 100.0 * n / total_called if total_called else 0
            print(f"  {base}: {n}  ({pct:.1f}% of base calls)")

    if args.window:
        start, end = (int(x) for x in args.window.split("-"))
        dist = net_length_change(args.bam, args.chrom, start, end)
        total = sum(dist.values())
        print(f"{args.chrom}:{start}-{end}  reads spanning window={total}")
        for net, n in dist.items():
            pct = 100.0 * n / total if total else 0
            label = "reference length" if net == 0 else f"{net:+d} bases"
            print(f"  {label}: {n}  ({pct:.1f}%)")


if __name__ == "__main__":
    main()
