#!/usr/bin/env python3
"""Check every genome, amplicon by amplicon, for stretches closer to a NUMT than to the
genome's own mitoLEAF haplogroup.

rtn judges one primer-trimmed read at a time, so the window is the insert of each NimaGen
amplicon (between its primers, located on rCRS from the primer files; the amplicon over the
origin is joined across 16569/1). Distance to the haplogroup = the genome's private and
missing differences (scan.tsv) inside the insert, a run of deleted or inserted bases counted
once. Distance to a NUMT = bwa mem edit distance to the nearest NUMT, unaligned ends counted
as differences. Only inserts at least 2 away from their haplogroup are checked; the others
cannot be clearly closer to a NUMT.

usage: amplicon_check.py outdir numts.fa resources_dir
       (outdir holds scan.tsv, aln.sam, singles.fa; resources_dir holds rCRS/ and primers/)
"""
import bisect
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

OUT, NUMTS, RES = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
TRUE_MT_LENGTH = 16569
IUPAC = {"R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT", "B": "CGT", "D": "AGT",
         "H": "ACT", "V": "ACG", "N": "ACGT"}


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


def load_primers(path):
    primers, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split("_")[0]
        elif line:
            primers.setdefault(name, []).append(line.strip("^$X").upper())
    return primers


def inserts():
    """{amplicon: [(start, end), ...]}: the rCRS segments between the primers (two for the
    amplicon over the origin)."""
    rcrs = "".join(l.strip() for l in open(RES / "rCRS/rCRS_NimaGen.fasta") if l[0] != ">")[:TRUE_MT_LENGTH].upper()
    circular = rcrs + rcrs[:500]
    regex = lambda p: re.compile("".join(f"[{IUPAC[c]}]" if c in IUPAC else c for c in p))
    left, right = load_primers(RES / "primers/left_primers.fasta"), load_primers(RES / "primers/right_primers_rc.fasta")
    out = {}
    for amp in left:
        lm = regex(left[amp][-1]).search(circular)
        rm = regex(right[amp][0]).search(circular, lm.end())
        start, end = lm.end() + 1, rm.start()  # 1-based, on the circular string
        out[amp] = [(start, end)] if end <= TRUE_MT_LENGTH else [(start, TRUE_MT_LENGTH), (1, end - TRUE_MT_LENGTH)]
    return out


def ref_pos(p):
    """rCRS position -> position on ref.fa (rCRS without the N at 3107)"""
    return p if p <= 3106 else p - 1


scan = pd.read_csv(OUT / "scan.tsv", sep="\t", keep_default_na=False)
ev = {r.genome: events(f"{r.private} {r.missing_list}".split()) for r in scan.itertuples()}
hg = dict(zip(scan.genome, scan.haplogroup))

# per genome: aligned blocks (ref start, query start, length), query on the original genome
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
            blocks[f[0]].append((r, q, n))
            q, r = q + n, r + n
        elif op in "IS":
            q += n
        elif op == "D":
            r += n
for b in blocks.values():
    b.sort()


def to_query(name, r, step):
    """Genome position of ref.fa position r; where r is deleted in the genome, the next
    aligned position in direction step."""
    b = blocks[name]
    for k in range(25):
        rr = r + k * step
        i = bisect.bisect_right(b, (rr, float("inf"), 0)) - 1
        if i >= 0 and b[i][0] <= rr < b[i][0] + b[i][2]:
            return b[i][1] + rr - b[i][0]
    return None


windows = inserts()
n_windows, checked = 0, {}
with open(OUT / "amplicons.fa", "w") as out:
    name = None
    for line in open(OUT / "singles.fa"):
        if line.startswith(">"):
            name = line[1:].strip()
            continue
        seq = line.strip()
        if name not in ev or name not in blocks:
            continue
        pos = ev[name]
        for amp, segs in sorted(windows.items()):
            n_windows += 1
            parts, d_hg = [], 0
            for lo, hi in segs:
                a, b = to_query(name, ref_pos(lo), +1), to_query(name, ref_pos(hi), -1)
                if a is None or b is None or not 0 <= b - a < 2 * (hi - lo + 1) + 50:
                    break  # not aligned there, or across a junction of glued pieces
                parts.append(seq[a:b + 1])
                d_hg += bisect.bisect_right(pos, hi) - bisect.bisect_left(pos, lo)
            else:
                if d_hg >= 2:
                    wid = f"{name}|{amp}"
                    checked[wid] = d_hg
                    out.write(f">{wid}\n{''.join(parts)}\n")
print("amplicon inserts", n_windows, "checked", len(checked))

sam = subprocess.run(["bwa", "mem", "-t", str(os.cpu_count()), NUMTS, str(OUT / "amplicons.fa")],
                     capture_output=True, text=True, check=True).stdout
d_numt, numt, length = {}, {}, {}
for line in sam.splitlines():
    if line[0] == "@":
        continue
    f = line.split("\t")
    if not int(f[1]) & 256:
        length.setdefault(f[0], len(f[9]))
    if int(f[1]) & (4 | 256):
        continue
    nm = int(re.search(r"NM:i:(\d+)", line)[1])
    clip = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in "SH")
    if f[0] not in d_numt or nm + clip < d_numt[f[0]]:
        d_numt[f[0]], numt[f[0]] = nm + clip, f"{f[2]}:{f[3]}"

rows = []
for wid, dh in checked.items():
    name, amp = wid.split("|")
    dn = d_numt.get(wid, length.get(wid, 0))  # not aligned to any NUMT: as far as its length
    rows.append({"genome": name, "haplogroup": hg[name], "amplicon": amp, "d_haplogroup": dh,
                 "d_numt": dn, "margin": dh - dn, "numt": numt.get(wid, "")})
t = pd.DataFrame(rows, columns=["genome", "haplogroup", "amplicon", "d_haplogroup", "d_numt", "margin", "numt"])
t.to_csv(OUT / "amplicons.tsv", sep="\t", index=False)
g = t.sort_values(["margin", "genome", "amplicon"], ascending=[False, True, True]).drop_duplicates("genome")
g.to_csv(OUT / "amplicon_genomes.tsv", sep="\t", index=False)
print("amplicon margins (d_haplogroup - d_numt):", t.margin.value_counts().sort_index().to_dict())
print("genomes by best margin:", g.margin.value_counts().sort_index().to_dict())
pd.set_option("display.width", 250)
print(g[g.margin >= 2].to_string(index=False))
