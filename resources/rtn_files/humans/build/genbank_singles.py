#!/usr/bin/env python3
"""Write the downloaded GenBank genomes one per line for scan_humans.py
(singles.fa, upper case, named by accession.version), and per genome its length,
the number of bases other than A, C, G, T and an md5 (sequence_stats.tsv).

usage: genbank_singles.py genbank.fa outdir
"""
import hashlib
import re
import sys
from pathlib import Path

OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
out = open(OUT / "singles.fa", "w")
stats = open(OUT / "sequence_stats.tsv", "w")
stats.write("genome\tlength\tnon_ACGT\tmd5\n")
name, parts = None, []


def flush():
    if name is None:
        return
    s = "".join(parts).upper()
    out.write(f">{name}\n{s}\n")
    stats.write(f"{name}\t{len(s)}\t{len(re.findall('[^ACGT]', s))}\t{hashlib.md5(s.encode()).hexdigest()}\n")


for line in open(sys.argv[1]):
    if line[0] == ">":
        flush()
        name, parts = line[1:].split()[0], []
    else:
        parts.append(line.strip())
flush()
