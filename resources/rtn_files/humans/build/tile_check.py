#!/usr/bin/env python3
"""Tile every humans.fa genome (200 bp, step 100) and flag tiles that are closer
to a NUMT than to the genome's own mitoLEAF haplogroup.

Distance to the haplogroup = the genome's private and missing differences
(scan.tsv) inside the tile's rCRS span, a run of deleted or inserted bases
counted once. Distance to a NUMT = bwa mem edit distance to the nearest NUMT,
unaligned tile ends counted as differences. Only tiles at least 2 away from
their haplogroup are checked; the others cannot be clearly closer to a NUMT.

usage: tile_check.py outdir numts.fa   (outdir holds scan.tsv, aln.sam, singles.fa)
"""
import bisect
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

OUT, NUMTS = Path(sys.argv[1]), sys.argv[2]
TILE, STEP = 200, 100
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def rcrs_pos(p):
    return p if p <= 3106 else p + 1


def events(tokens):
    """rCRS positions of difference events; consecutive deleted bases and the
    bases of one insertion count once."""
    out, last_del, last_ins = [], None, None
    for t in tokens:
        p = int(re.match(r"\d+", t)[0])
        if t.endswith("-"):
            if last_del is None or p != last_del + 1:
                out.append(p)
            last_del = p
        elif "." in t:
            if p != last_ins:
                out.append(p)
            last_ins = p
        else:
            out.append(p)
    return sorted(out)


scan = pd.read_csv(OUT / "scan.tsv", sep="\t", keep_default_na=False)
ev = {r.genome: events(f"{r.private} {r.missing_list}".split()) for r in scan.itertuples()}
hg = dict(zip(scan.genome, scan.haplogroup))

# per genome: aligned blocks (query start, rCRS-without-N start, length), query on the original genome
blocks = defaultdict(list)
for line in open(OUT / "aln.sam"):
    if line[0] == "@":
        continue
    f = line.split("\t")
    flag = int(f[1])
    if flag & (4 | 256) or flag & 16:
        continue
    ops = [(int(n), op) for n, op in CIGAR.findall(f[5])]
    q = ops[0][0] if ops[0][1] == "H" else 0  # hard-clipped lead of supplementary pieces
    r = int(f[3])
    for n, op in ops:
        if op == "M":
            blocks[f[0]].append((q, r, n))
            q, r = q + n, r + n
        elif op in "IS":
            q += n
        elif op == "D":
            r += n
for b in blocks.values():
    b.sort()


def to_ref(name, q):
    b = blocks[name]
    i = bisect.bisect_right(b, (q, float("inf"), 0)) - 1
    if i >= 0 and b[i][0] <= q < b[i][0] + b[i][2]:
        return b[i][1] + q - b[i][0]
    return None


n_tiles, tiles = 0, {}
with open(OUT / "tiles.fa", "w") as out:
    name = None
    for line in open(OUT / "singles.fa"):
        if line.startswith(">"):
            name = line[1:].strip()
            continue
        seq = line.strip()
        starts = list(range(0, max(len(seq) - TILE, 0) + 1, STEP))
        if starts and starts[-1] + TILE < len(seq):
            starts.append(len(seq) - TILE)
        pos = ev[name]
        for s in starts:
            n_tiles += 1
            a = next((to_ref(name, q) for q in range(s, s + 20) if to_ref(name, q)), None)
            b = next((to_ref(name, q) for q in range(s + TILE - 1, s + TILE - 21, -1) if to_ref(name, q)), None)
            if a is None or b is None or not 0 < b - a < 2 * TILE:
                continue  # unaligned or across a junction of glued pieces
            lo, hi = rcrs_pos(a), rcrs_pos(b)
            d_hg = bisect.bisect_right(pos, hi) - bisect.bisect_left(pos, lo)
            if d_hg >= 2:
                tid = f"{name}|{s}|{lo}-{hi}"
                tiles[tid] = d_hg
                out.write(f">{tid}\n{seq[s:s + TILE]}\n")
print("tiles", n_tiles, "checked", len(tiles))

sam = subprocess.run(["bwa", "mem", "-t", str(os.cpu_count()), NUMTS, str(OUT / "tiles.fa")],
                     capture_output=True, text=True, check=True).stdout
d_numt, numt = {}, {}
for line in sam.splitlines():
    if line[0] == "@":
        continue
    f = line.split("\t")
    if int(f[1]) & (4 | 256):
        continue
    nm = int(re.search(r"NM:i:(\d+)", line)[1])
    clip = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in "SH")
    if f[0] not in d_numt or nm + clip < d_numt[f[0]]:
        d_numt[f[0]], numt[f[0]] = nm + clip, f"{f[2]}:{f[3]}"

rows = []
for tid, dh in tiles.items():
    dn = d_numt.get(tid, TILE)
    name, s, span = tid.split("|")
    rows.append({"genome": name, "haplogroup": hg[name], "tile_start": int(s), "rCRS_span": span,
                 "d_haplogroup": dh, "d_numt": dn, "margin": dh - dn, "numt": numt.get(tid, "")})
t = pd.DataFrame(rows)
t.to_csv(OUT / "tiles.tsv", sep="\t", index=False)
g = t.sort_values("margin", ascending=False).drop_duplicates("genome")
g.to_csv(OUT / "tile_genomes.tsv", sep="\t", index=False)
print("tile margins (d_haplogroup - d_numt):", t.margin.value_counts().sort_index().to_dict())
print("genomes by best margin:", g.margin.value_counts().sort_index().to_dict())
pd.set_option("display.width", 250)
print(g[g.margin >= 2].to_string(index=False))
