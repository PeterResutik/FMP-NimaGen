#!/usr/bin/env python3
"""For outlier genomes: cut the genome around its densest private cluster
(rCRS window +-100 bp, mapped through the genome's alignment) and compare
its distance to rCRS with its distance to the nearest NUMT.

usage: numt_check.py outdir numts.fa   (outdir holds scan.tsv, aln.sam, singles.fa, ref.fa)
"""
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

OUT, NUMTS = Path(sys.argv[1]), sys.argv[2]
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def no_n(p):
    """rCRS position -> position on rCRS without its N at 3107."""
    return p if p <= 3106 else p - 1


d = pd.read_csv(OUT / "scan.tsv", sep="\t")
d["pieces"] = d.covered.str.count(";") + 1
cand = d[((d.cluster50 >= 4) | (d.private_subs >= 25)) & (d.pieces < 20)].copy()
names = set(cand.genome)

seqs, name = {}, None
for line in open(OUT / "singles.fa"):
    if line.startswith(">"):
        name = line[1:].strip()
    elif name in names:
        seqs[name] = line.strip()

# rCRS (no N) position -> query position, from each genome's alignments
qmap = {n: {} for n in names}
for line in open(OUT / "aln.sam"):
    if line[0] == "@":
        continue
    f = line.split("\t")
    if f[0] not in names or int(f[1]) & (4 | 256):
        continue
    r, q = int(f[3]), 0
    ops = CIGAR.findall(f[5])
    if int(f[1]) & 2048:
        continue  # supplementary pieces are hard-clipped; the primary covers the clusters here
    for n, op in ops:
        n = int(n)
        if op in "M":
            for k in range(n):
                qmap[f[0]].setdefault(r + k, q + k)
            r, q = r + n, q + n
        elif op in "IS":
            q += n
        elif op == "D":
            r += n

with open(OUT / "segments.fa", "w") as out:
    for _, row in cand.iterrows():
        pos = [int(re.match(r"\d+", t)[0]) for t in str(row.cluster).split()]
        if not pos:
            continue
        a, b = no_n(max(1, pos[0] - 100)), no_n(min(16569, pos[-1] + 100))
        m = qmap[row.genome]
        qa = next((m[p] for p in range(a, b) if p in m), None)
        qb = next((m[p] for p in range(b, a, -1) if p in m), None)
        if qa is None or qb is None or qb <= qa:
            continue
        out.write(f">{row.genome}\n{seqs[row.genome][qa:qb + 1]}\n")


def best_nm(index):
    sam = subprocess.run(["bwa", "mem", "-a", str(index), str(OUT / "segments.fa")],
                         capture_output=True, text=True, check=True).stdout
    best = {}
    for line in sam.splitlines():
        if line[0] == "@":
            continue
        f = line.split("\t")
        if int(f[1]) & 4:
            continue
        nm = int(re.search(r"NM:i:(\d+)", line)[1])
        clip = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in "SH")
        score = nm + clip  # unaligned ends count as differences
        if f[0] not in best or score < best[f[0]][0]:
            best[f[0]] = (score, f[2])
    return best


if not Path(NUMTS + ".bwt").exists():
    subprocess.run(["bwa", "index", NUMTS], check=True, capture_output=True)
to_rcrs, to_numt = best_nm(OUT / "ref.fa"), best_nm(NUMTS)
cand["segment_vs_rCRS"] = cand.genome.map(lambda g: to_rcrs.get(g, (None,))[0])
cand["segment_vs_NUMT"] = cand.genome.map(lambda g: to_numt.get(g, (None,))[0])
cand["nearest_NUMT"] = cand.genome.map(lambda g: to_numt.get(g, (None, None))[1])
cand["closer_to_NUMT"] = cand.segment_vs_NUMT < cand.segment_vs_rCRS
cols = ["genome", "length", "haplogroup", "private_subs", "private_indel_bases", "missing", "cluster50",
        "segment_vs_rCRS", "segment_vs_NUMT", "nearest_NUMT", "closer_to_NUMT", "cluster"]
cand[cols].sort_values(["closer_to_NUMT", "cluster50", "private_subs"], ascending=False).to_csv(
    OUT / "numt_check.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 70)
print(cand[cols].sort_values(["closer_to_NUMT", "cluster50", "private_subs"], ascending=False).to_string(index=False))
