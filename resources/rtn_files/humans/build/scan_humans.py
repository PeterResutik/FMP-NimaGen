#!/usr/bin/env python3
"""Outlier scan of rtn's humans.fa against the mitoLEAF tree.

Each genome (one copy) is aligned to rCRS, its differences are written in
mitoLEAF's notation (substitutions 73G, deletions per base 3'-placed 8281-,
insertions 573.1C), and it is matched to the haplogroup whose motif explains
most of them. Per genome: the differences the motif does not explain
(private), the motif variants the genome lacks (missing), and the densest
50-bp cluster of private substitutions.

usage: scan_humans.py humans.fa hgmotifs.json outdir
"""
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HUMANS, MOTIFS, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
# indels here are length variation, not counted as private or missing
LENGTH_REGIONS = [(57, 60), (300, 316), (452, 463), (514, 525), (568, 573), (956, 965),
                  (5895, 5899), (8272, 8289), (16180, 16193), (3105, 3110)]
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
REF = None


def rcrs_pos(p):
    """Position on rCRS without its N at 3107 -> rCRS numbering."""
    return p if p <= 3106 else p + 1


def tpos(token):
    return int(re.match(r"\d+", token)[0])


def in_length_region(token):
    p = tpos(token)
    return ("-" in token or "." in token) and any(a <= p <= b for a, b in LENGTH_REGIONS)


def split_humans():
    singles, ref = OUT / "singles.fa", OUT / "ref.fa"
    if singles.exists():
        return
    with open(HUMANS) as f, open(singles, "w") as out:
        name, parts = None, []

        def flush():
            s = "".join(parts)
            h = len(s) // 2
            out.write(f">{name}\n{s[:h]}\n")
            if name == "chrM":
                ref.write_text(f">chrM\n{s[:h]}\n")

        for line in f:
            if line.startswith(">"):
                if name:
                    flush()
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
        flush()


def align():
    sam = OUT / "aln.sam"
    if sam.exists():
        return
    subprocess.run(["bwa", "index", str(OUT / "ref.fa")], check=True, capture_output=True)
    with open(sam, "w") as f:
        subprocess.run(["bwa", "mem", "-t", str(os.cpu_count()), str(OUT / "ref.fa"), str(OUT / "singles.fa")],
                       stdout=f, stderr=subprocess.DEVNULL, check=True)


def init(ref):
    global REF
    REF = ref


def record_tokens(line):
    f = line.split("\t")
    flag = int(f[1])
    if flag & 4 or flag & 256:
        return f[0], [], None
    start, seq = int(f[3]), f[9]
    r, q, toks = start - 1, 0, []
    for n, op in CIGAR.findall(f[5]):
        n = int(n)
        if op == "M":
            a = np.frombuffer(REF[r:r + n].encode(), dtype=np.uint8)
            b = np.frombuffer(seq[q:q + n].encode(), dtype=np.uint8)
            for i in np.nonzero(a != b)[0]:
                toks.append(f"{rcrs_pos(r + i + 1)}{seq[q + i]}")
            r, q = r + n, q + n
        elif op == "I":
            p, s = r, seq[q:q + n]
            while p < len(REF) and REF[p] == s[0]:  # shift 3': REF[p] is position p+1
                s, p = s[1:] + s[0], p + 1
            toks += [f"{rcrs_pos(p)}.{k + 1}{c}" for k, c in enumerate(s)]
            q += n
        elif op == "D":
            p, e = r + 1, r + n
            while e < len(REF) and REF[e] == REF[p - 1]:  # shift 3'
                p, e = p + 1, e + 1
            toks += [f"{rcrs_pos(x)}-" for x in range(p, e + 1)]
            r += n
        elif op == "S":
            q += n
    return f[0], toks, (rcrs_pos(start), rcrs_pos(r))


def densest(positions, width=50):
    best, window = 0, (0, 0)
    j = 0
    for i in range(len(positions)):
        while positions[i] - positions[j] >= width:
            j += 1
        if i - j + 1 > best:
            best, window = i - j + 1, (j, i + 1)
    return best, window


def main():
    split_humans()
    align()
    ref = "".join(l.strip() for l in open(OUT / "ref.fa") if not l.startswith(">"))
    lengths = {}
    for line in open(OUT / "singles.fa"):
        if line.startswith(">"):
            name = line[1:].strip()
        else:
            lengths[name] = len(line.strip())
    tokens, cover = defaultdict(set), defaultdict(list)
    with open(OUT / "aln.sam") as f, Pool(os.cpu_count(), initializer=init, initargs=(ref,)) as pool:
        lines = [l for l in f if not l.startswith("@")]
        for name, toks, iv in pool.imap(record_tokens, lines, chunksize=200):
            tokens[name].update(toks)
            if iv:
                cover[name].append(iv)

    motifs = json.load(open(MOTIFS))
    hgs = list(motifs)
    fixed, ambiguous = [], []
    for hg in hgs:
        fx, amb = set(), set()
        for t in motifs[hg].split():
            suffix = re.fullmatch(r"\d+(?:\.\d+)?([A-Za-z-])", t)[1]
            (fx.add(t) if suffix in "ACGT-" else amb.add(tpos(t)))
        fixed.append(fx)
        ambiguous.append(amb)
    vocab = sorted(set().union(*fixed), key=lambda t: (tpos(t), t))
    index = {t: i for i, t in enumerate(vocab)}
    vpos = np.array([tpos(t) for t in vocab])
    M = np.zeros((len(hgs), len(vocab)), dtype=np.float32)
    for h, fx in enumerate(fixed):
        M[h, [index[t] for t in fx]] = 1

    genomes = list(lengths)
    G = np.zeros((len(genomes), len(vocab)), dtype=np.float32)
    C = np.zeros((len(genomes), len(vocab)), dtype=np.float32)
    for g, name in enumerate(genomes):
        hits = [index[t] for t in tokens[name] if t in index]
        G[g, hits] = 1
        for a, b in cover[name]:
            C[g, (vpos >= a) & (vpos <= b)] = 1
    best = np.zeros(len(genomes), dtype=int)
    for s in range(0, len(genomes), 4000):
        inter = G[s:s + 4000] @ M.T
        mcov = C[s:s + 4000] @ M.T
        best[s:s + 4000] = np.argmax(2 * inter - mcov, axis=1)

    rows = []
    for g, name in enumerate(genomes):
        h = best[g]
        fx, amb = fixed[h], ambiguous[h]
        own = {t for t in tokens[name] if tpos(t) not in amb}
        private = sorted((t for t in own - fx if not in_length_region(t)), key=lambda t: (tpos(t), t))
        covered = lambda p: any(a <= p <= b for a, b in cover[name])
        missing = sorted((t for t in fx - own if covered(tpos(t)) and not in_length_region(t)),
                         key=lambda t: (tpos(t), t))
        subs = [t for t in private if re.fullmatch(r"\d+[ACGT]", t)]
        indels = [t for t in private if t not in subs]
        n, (i, j) = densest([tpos(t) for t in subs])
        rows.append({"genome": name, "length": lengths[name],
                     "covered": ";".join(f"{a}-{b}" for a, b in sorted(cover[name])),
                     "haplogroup": hgs[h], "differences": len(tokens[name]),
                     "private_subs": len(subs), "private_indel_bases": len(indels), "missing": len(missing),
                     "cluster50": n, "cluster": " ".join(subs[i:j]),
                     "private": " ".join(private), "missing_list": " ".join(missing)})
    rows.sort(key=lambda r: (-r["cluster50"], -r["private_subs"]))
    cols = list(rows[0])
    with open(OUT / "scan.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("genomes", len(rows), "unaligned", sum(1 for r in rows if not r["covered"]))


if __name__ == "__main__":
    main()
