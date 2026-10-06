#!/usr/bin/env python3
"""Rank mtDNA positions by how often they change in the human phylogeny.

First the hypervariable positions 16519, 16182 and 16183, which PhyloTree and
mitoLEAF leave out of the tree by convention; then every other position by the
number of branches of mitoLEAF's tree with a substitution there, forward or back
(insertions, deletions and ambiguity codes not counted; ties by position).

usage: rank_positions.py mitoLEAF_phm.xml N > recurrent_positions.tsv
"""
import re
import sys
from collections import Counter

HYPERVARIABLE = [16519, 16182, 16183]


def branch_counts(xml):
    """{position: number of branches with a substitution there}; the root's motif is not a branch"""
    counts = Counter()
    for m in re.finditer(r'<Node Id="([^"]+)" HG="([^"]*)"', xml):
        if m[1] == "mt-MRCA":
            continue
        for t in m[2].split():
            s = re.fullmatch(r"(\d+)[ACGT]", t) or re.fullmatch(r"[ACGT](\d+)", t)
            if s:
                counts[int(s[1])] += 1
    return counts


def ranking(counts):
    """positions in rank order: HYPERVARIABLE first, then by count (most first) and position"""
    rest = sorted((p for p in counts if p not in HYPERVARIABLE), key=lambda p: (-counts[p], p))
    return HYPERVARIABLE + rest


if __name__ == "__main__":
    counts = branch_counts(open(sys.argv[1]).read())
    print("rank\tposition\tbranches")
    for i, p in enumerate(ranking(counts)[:int(sys.argv[2])], 1):
        print(f"{i}\t{p}\t{counts.get(p, 0)}")
