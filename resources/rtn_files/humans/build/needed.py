#!/usr/bin/env python3
"""Which mitoLEAF haplogroups add something to humans.fa for NimaGen reads?

Per haplogroup and amplicon (library range, primers included): differences
between the haplogroup's motif and the nearest complete humans.fa genome over
that amplicon (positions ambiguous in most motifs and indels in length regions
left out; the excluded genomes left out). Writes needed.tsv; the haplogroups
whose worst amplicon is 1 or more differences away are the ones worth adding.

usage: needed.py outdir hgmotifs.json library.txt exclude.txt   (outdir holds scan.tsv, aln.sam, ref.fa)
"""
import json
import os
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

OUT, MOTIFS, LIBRARY = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
EXCLUDED = {l.split()[0] for l in open(sys.argv[4]) if l.strip() and not l.startswith("#")}
PRIMER = 30  # bases added on each side of a library range for the primers
LENGTH_REGIONS = [(57, 60), (300, 316), (452, 463), (514, 525), (568, 573), (956, 965),
                  (5895, 5899), (8272, 8289), (16180, 16193), (3105, 3110)]
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF = None


def rcrs_pos(p):
    return p if p <= 3106 else p + 1


def tpos(t):
    return int(re.match(r"\d+", t)[0])


def init(ref):
    global REF
    REF = ref


def record_tokens(line):
    f = line.split("\t")
    if int(f[1]) & (4 | 256):
        return f[0], []
    r, q, seq, toks = int(f[3]) - 1, 0, f[9], []
    for n, op in CIGAR.findall(f[5]):
        n = int(n)
        if op == "M":
            toks += [f"{rcrs_pos(r + i + 1)}{seq[q + i]}" for i in range(n) if REF[r + i] != seq[q + i]]
            r, q = r + n, q + n
        elif op == "I":
            p, s = r, seq[q:q + n]
            while p < len(REF) and REF[p] == s[0]:
                s, p = s[1:] + s[0], p + 1
            toks += [f"{rcrs_pos(p)}.{k + 1}{c}" for k, c in enumerate(s)]
            q += n
        elif op == "D":
            p, e = r + 1, r + n
            while e < len(REF) and REF[e] == REF[p - 1]:
                p, e = p + 1, e + 1
            toks += [f"{rcrs_pos(x)}-" for x in range(p, e + 1)]
            r += n
        elif op == "S":
            q += n
    return f[0], toks


def windows():
    out = {}
    for line in open(LIBRARY):
        m = re.match(r"(mtNG_\d+) = chrM, (\d+), (\d+)(?:, (\d+), (\d+))?", line)
        if m and m[1] not in out:
            if m[4]:  # origin amplicon: two ranges
                out[m[1]] = [(int(m[2]) - PRIMER, 16569), (1, int(m[5]) + PRIMER)]
            else:
                out[m[1]] = [(max(1, int(m[2]) - PRIMER), min(16569, int(m[3]) + PRIMER))]
    return out


def main():
    scan = pd.read_csv(OUT / "scan.tsv", sep="\t", keep_default_na=False)
    complete = set(scan.genome[(scan.covered.str.count(";") == 0) & (scan.length >= 16500)]) - EXCLUDED
    ref = "".join(l.strip() for l in open(OUT / "ref.fa") if not l.startswith(">"))
    tokens = defaultdict(set)
    with open(OUT / "aln.sam") as f, Pool(os.cpu_count(), initializer=init, initargs=(ref,)) as pool:
        lines = [l for l in f if l[0] != "@" and l.split("\t", 1)[0] in complete]
        for name, toks in pool.imap(record_tokens, lines, chunksize=200):
            tokens[name].update(toks)

    motifs = json.load(open(MOTIFS))
    hgs = list(motifs)
    amb = Counter(tpos(t) for m in motifs.values() for t in m.split()
                  if not re.fullmatch(r"\d+(?:\.\d+)?[ACGT-]", t))
    common_amb = {p for p, c in amb.items() if c > len(hgs) / 2}
    keep = lambda t: tpos(t) not in common_amb and not (
        ("-" in t or "." in t) and any(a <= tpos(t) <= b for a, b in LENGTH_REGIONS))
    fixed = [{t for t in motifs[h].split() if re.fullmatch(r"\d+(?:\.\d+)?[ACGT-]", t) and keep(t)} for h in hgs]
    genomes = sorted(complete)
    gtok = [{t for t in tokens[g] if keep(t)} for g in genomes]

    wins = windows()
    mind = np.zeros((len(hgs), len(wins)), dtype=np.int32)
    for w, (amp, ranges) in enumerate(sorted(wins.items())):
        inside = lambda t: any(a <= tpos(t) <= b for a, b in ranges)
        hw = [{t for t in fx if inside(t)} for fx in fixed]
        gw = [{t for t in g if inside(t)} for g in gtok]
        vocab = sorted(set().union(*hw))
        index = {t: i for i, t in enumerate(vocab)}
        H = np.zeros((len(hgs), max(len(vocab), 1)), dtype=np.float32)
        for h, s in enumerate(hw):
            H[h, [index[t] for t in s]] = 1
        G = np.zeros((len(genomes), max(len(vocab), 1)), dtype=np.float32)
        for g, s in enumerate(gw):
            G[g, [index[t] for t in s if t in index]] = 1
        gsize = np.array([len(s) for s in gw], dtype=np.float32)
        best = np.full(len(hgs), np.inf, dtype=np.float32)
        for s in range(0, len(genomes), 8000):
            d = H.sum(1)[:, None] + gsize[None, s:s + 8000] - 2 * (H @ G[s:s + 8000].T)
            best = np.minimum(best, d.min(1))
        mind[:, w] = best
    amps = sorted(wins)
    worst = mind.max(1)
    far = (mind >= 2).sum(1)
    rows = pd.DataFrame({"haplogroup": hgs, "worst_amplicon_distance": worst, "amplicons_2plus": far,
                         "amplicons": [" ".join(f"{amps[w]}:{mind[h, w]}" for w in np.nonzero(mind[h] >= 2)[0])
                                       for h in range(len(hgs))]})
    rows.to_csv(OUT / "needed.tsv", sep="\t", index=False)
    print("haplogroups:", len(hgs), " amplicons:", len(amps))
    print("worst amplicon distance to the nearest genome:",
          pd.Series(worst).clip(upper=6).value_counts().sort_index().to_dict(), "(6 = 6 or more)")
    print("needed (some amplicon 2+ from every genome):", int((worst >= 2).sum()),
          " of which some amplicon 3+:", int((worst >= 3).sum()))
    per_amp = (mind >= 2).sum(0)
    print("amplicons most often uncovered:", sorted(zip(per_amp, amps), reverse=True)[:10])
    for hg in ["L7a1a1a", "L3f3", "L0d2b1b1", "L2b1a3a", "L3b1a12", "L1c3a1a", "H2a2a1", "U5a1a1"]:
        if hg in hgs:
            r = rows[rows.haplogroup == hg].iloc[0]
            print(f"  {hg}: worst {r.worst_amplicon_distance}, amplicons 2+: {r.amplicons}")


if __name__ == "__main__":
    main()
