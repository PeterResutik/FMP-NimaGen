#!/usr/bin/env python3
"""Base qualities for rtn: lower them at chosen positions in the copy of the reads rtn
evaluates, and give rtn's output the reads' own qualities back.

usage: samtools view -h in.bam | rtn_qualities.py lower 152,16311,... > lowq.sam
       samtools view -h lowq.rtn.bam | rtn_qualities.py restore in.bam > restored.sam

lower:   base quality 2 ('#') at every read base aligned to one of the positions (rCRS
         numbering; the appended copy of 1-53, at 16570-16622, counts as 1-53).
restore: each record takes the qualities of the record in in.bam with the same name, the
         same secondary/supplementary flags and the same position; rtn's MAPQ and tags stay.
         A record without its original qualities is an error.
"""
import re
import subprocess
import sys

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')
LOW = '#'
RCRS_LENGTH = 16569


def lower(line, positions):
    """one SAM line with base quality 2 at the read bases aligned to positions"""
    if line[0] == '@':
        return line
    f = line.split('\t')
    if int(f[1]) & 4 or f[10] == '*':
        return line
    qual = list(f[10])
    r, q = int(f[3]), 0
    for n, op in CIGAR.findall(f[5]):
        n = int(n)
        if op in 'M=X':
            for i in range(n):
                p = r + i
                if (p if p <= RCRS_LENGTH else p - RCRS_LENGTH) in positions:
                    qual[q + i] = LOW
            r, q = r + n, q + n
        elif op in 'IS':
            q += n
        elif op in 'DN':
            r += n
    f[10] = ''.join(qual)
    return '\t'.join(f)


def key(f):
    return f[0], int(f[1]) & 0x900, f[3]


def originals(bam):
    """{(name, secondary/supplementary flags, position): qualities} of the records in bam"""
    out = {}
    view = subprocess.run(['samtools', 'view', bam], capture_output=True, text=True, check=True)
    for line in view.stdout.splitlines():
        f = line.split('\t', 11)
        out[key(f)] = f[10]
    return out


def restore(line, orig):
    """one SAM line with its original qualities, or None if there are none of the same length"""
    if line[0] == '@':
        return line
    f = line.rstrip('\n').split('\t')
    q = orig.get(key(f))
    if q is None or len(q) != len(f[10]):
        return None
    f[10] = q
    return '\t'.join(f) + '\n'


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ('lower', 'restore'):
        sys.exit(__doc__)
    if sys.argv[1] == 'lower':
        positions = {int(p) for p in sys.argv[2].split(',')}
        for line in sys.stdin:
            sys.stdout.write(lower(line, positions))
        return
    orig, missing = originals(sys.argv[2]), 0
    for line in sys.stdin:
        out = restore(line, orig)
        if out is None:
            missing += 1
        else:
            sys.stdout.write(out)
    if missing:
        sys.exit(f'rtn_qualities.py restore: {missing} records without their original qualities')


if __name__ == '__main__':
    main()
