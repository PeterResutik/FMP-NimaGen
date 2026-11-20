#!/usr/bin/env python3
"""Every mapped read of a NUMT-filter step in one BAM: the reads the step removes get MAPQ 0, and every record gets
a tag ZF with the reason, so that the removed reads can be seen next to the kept ones (in IGV: colour or group the
alignments by the tag ZF). Unmapped reads are left out (IGV cannot show them).

usage: samtools view -h in.bam | mark_removed_reads.py p09 in.rtn.bam MAPQ CLIPPED removed.tsv > all_reads.sam
       samtools view -h merged.bam | mark_removed_reads.py p10 kept_names.txt removed.tsv MAPQ > all_reads.sam

ZF values:
  kept      the read is used
  numt      rtn removed the read (set its MAPQ to 0): no human mitogenome explains it closely enough
  clipped   the mapping without primers could not align the read end to end (only when CLIPPED is 1)
  mapq      bwa gave the read a mapping quality below MAPQ, and rtn did not remove it (rtn does not judge reads
            below its own minimum, MAPQ 4)
  trimmed   (p10) primer trimming removed the read, so rtn never judged it
  unmapped  (p10) the read maps with its primers but not without them, so rtn never judged it

p09: every mapped record of in.bam (the reads without primers), as rtn wrote it (with its tags ZH and ZN) where
     rtn wrote it; removed.tsv gets the reason of each removed read, from its primary record (unmapped reads too).
p10: every mapped record of the reads with primers; a read p09 kept is 'kept', or 'mapq' when its own mapping quality
     is below MAPQ; any other read takes its reason from removed.tsv, or 'trimmed' if p09 never had it.
"""
import subprocess
import sys

TAG = 'ZF:Z:'


def key(f):
    return f[0], int(f[1]) & 0x900, f[3]


def is_primary(f):
    return not int(f[1]) & 0x900


def is_unmapped(f):
    return bool(int(f[1]) & 4)


def is_clipped(f):
    return 'S' in f[5] or 'H' in f[5]


def mark(f, reason):
    """one SAM line: MAPQ 0 unless kept, and the tag ZF"""
    f = list(f)
    if reason != 'kept':
        f[4] = '0'
    return '\t'.join(f + [TAG + reason]) + '\n'


def p09_reason(f, rtn, mapq, clipped):
    """the reason for one record of the reads without primers, given rtn's record of it (None: rtn left it out)"""
    if is_unmapped(f):
        return 'unmapped'
    if rtn is not None and rtn[4] != f[4]:  # rtn sets the reads it removes to MAPQ 0
        return 'numt'
    if int(f[4]) < mapq:
        return 'mapq'
    if rtn is None:
        return None
    if clipped and is_clipped(f):
        return 'clipped'
    return 'kept'


def p09(lines, rtn_records, mapq, clipped):
    """(SAM lines of the mapped reads, {name: reason} of the removed reads)"""
    out, records = [], []
    for line in lines:
        if line[0] == '@':
            out.append(line)
            continue
        f = line.rstrip('\n').split('\t')
        rtn = rtn_records.get(key(f))
        records.append((f, rtn, p09_reason(f, rtn, mapq, clipped)))
    primary = {f[0]: reason for f, _, reason in records if is_primary(f)}
    removed = {}
    for f, rtn, reason in records:
        if reason is None:  # a mapped record rtn left out (e.g. the supplementary part of a split read)
            reason = primary.get(f[0])
            reason = 'numt' if reason in (None, 'kept') else reason
        if reason != 'unmapped':
            out.append(mark(rtn or f, reason))
        if is_primary(f) and reason != 'kept':
            removed[f[0]] = reason
    return out, removed


def p10_reason(f, kept, removed, mapq):
    """the reason for one mapped record of the reads with primers"""
    if f[0] in kept:
        return 'mapq' if int(f[4]) < mapq else 'kept'
    return removed.get(f[0], 'trimmed')


def records(bam):
    """{(name, secondary/supplementary flags, position): fields} of the records in bam"""
    view = subprocess.run(['samtools', 'view', bam], capture_output=True, text=True, check=True)
    return {key(f): f for f in (line.split('\t') for line in view.stdout.splitlines())}


def main():
    a = sys.argv[1:]
    if a[:1] == ['p09'] and len(a) == 5:
        out, removed = p09(sys.stdin, records(a[1]), int(a[2]), a[3] == '1')
        sys.stdout.writelines(out)
        with open(a[4], 'w') as h:
            h.writelines(f'{name}\t{reason}\n' for name, reason in sorted(removed.items()))
    elif a[:1] == ['p10'] and len(a) == 4:
        kept = {line.strip() for line in open(a[1]) if line.strip()}
        removed = dict(line.rstrip('\n').split('\t') for line in open(a[2]) if line.strip())
        mapq = int(a[3])
        for line in sys.stdin:
            if line[0] == '@':
                sys.stdout.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if not is_unmapped(f):
                sys.stdout.write(mark(f, p10_reason(f, kept, removed, mapq)))
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
