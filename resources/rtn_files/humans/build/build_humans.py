#!/usr/bin/env python3
"""Assemble the human mitogenomes for rtn from the screened sources.

Order of the entries: chrM (rCRS), every mitoLEAF haplogroup, the GenBank genomes of
mitoLEAF's list that pass the screen, the re-called 1000 Genomes and HGDP people that
pass it. A sequence that occurs more than once is written once, under the first name;
the other names go to members.tsv. Every entry is followed by its own first 100 bases,
so reads across the origin align in one piece.

A genome is left out (left_out.tsv, with the reason) when:
- it is not aligned to rCRS in one piece;
- the insert of a NimaGen amplicon (what rtn judges in a primer-trimmed read) is closer to
  a NUMT than to the genome's own haplogroup by 3 or more (amplicon_check.py);
- it has more than 10 unconfirmed substitutions, or 4 or more of them within 50 bases.
  A substitution is confirmed when a genome outside the genome's own accession series
  (same prefix, number in the same block of 500) carries it too; the genomes counted
  are those of mitoLEAF's list and the re-called people;
- (GenBank) it has more than 5 bases other than A, C, G, T;
- (GenBank) the accession belongs to a submission listed in SERIES.
One to five bases other than A, C, G, T are written as the rCRS base at that position.
The people's insertions and deletions in the length regions are not used (whole-genome
calls are unreliable there), so those stretches have the rCRS length.

usage: build_humans.py workdir people_calls.tsv.gz out.fa
  workdir holds ref.fa, mitoleaf.fa, scan/ and scan_people/ (scan_humans.py, amplicon_check.py)
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

import people

W, CALLS, OUT = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
TAIL = 100
SERIES = {"MN6281": "submission MN628108-MN628167 (Zhang 2019): 27 of its 69 genomes carry a stretch "
                    "identical to a NUMT on chromosome 1 or 3"}
UNCONFIRMED_MAX, WINDOW, WINDOW_MAX = 10, 50, 4
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
MD = re.compile(r"(\d+)|(\^[A-Z]+)|([A-Z])")
ref = "".join(l.strip() for l in open(W / "ref.fa") if l[0] != ">")


def substitutions(sam):
    """{genome: set of (position, base)} from the alignments to rCRS (MD tag)"""
    out = {}
    for line in open(sam):
        if line[0] == "@":
            continue
        f = line.split("\t")
        if int(f[1]) & 0x904:
            continue
        md = next((x[5:] for x in f[11:] if x.startswith("MD:Z:")), None)
        if md is None:
            continue
        qpos, r, q = {}, int(f[3]), 0
        for n, op in CIGAR.findall(f[5]):
            n = int(n)
            if op == "M":
                for i in range(n):
                    qpos[r + i] = q + i
                r, q = r + n, q + n
            elif op in "IS":
                q += n
            elif op == "D":
                r += n
        r, subs = int(f[3]), set()
        for num, dele, mm in MD.findall(md.strip()):
            if num:
                r += int(num)
            elif dele:
                r += len(dele) - 1
            else:
                b = f[9][qpos[r]] if r in qpos else "N"
                if b in "ACGT":
                    subs.add((r, b))
                r += 1
        out.setdefault(f[0], set()).update(subs)
    return out


def series(g):
    m = re.match(r"([A-Z]+)(\d+)", g)
    return (m.group(1), int(m.group(2)) // 500) if m else (g, 0)


SUBS = {**substitutions(W / "scan/aln.sam"), **substitutions(W / "scan_people/aln.sam")}
CARRIERS = defaultdict(set)
for g, v in SUBS.items():
    for x in v:
        CARRIERS[x].add(series(g))


def unconfirmed(g):
    """substitutions of g that no genome outside its own series carries; the most of them within WINDOW bases"""
    u = sorted(p for p, b in SUBS.get(g, ()) if not (CARRIERS[(p, b)] - {series(g)}))
    best, j = 0, 0
    for i in range(len(u)):
        while u[i] - u[j] > WINDOW:
            j += 1
        best = max(best, i - j + 1)
    return u, best


def screen(scan_dir, stats=None):
    s = pd.read_csv(scan_dir / "scan.tsv", sep="\t", keep_default_na=False).set_index("genome")
    a = pd.read_csv(scan_dir / "amplicon_genomes.tsv", sep="\t", keep_default_na=False,
                    dtype={"amplicon": str}).set_index("genome")
    m = a.margin.reindex(s.index).fillna(-99)
    why = {}
    for g in s.index:
        r = []
        if s.covered[g].count(";") > 0:
            r.append("not aligned in one piece")
        if m[g] >= 3:
            r.append(f"amplicon {a.amplicon[g]} closer to a NUMT than to its haplogroup by {int(m[g])}")
        u, w = unconfirmed(g)
        if len(u) > UNCONFIRMED_MAX:
            r.append(f"{len(u)} substitutions no genome outside its series carries")
        if w >= WINDOW_MAX:
            r.append(f"{w} such substitutions within {WINDOW} bases")
        if stats is not None and stats.non_ACGT[g] > 5:
            r.append(f"{stats.non_ACGT[g]} bases other than A, C, G, T")
        for prefix, reason in SERIES.items():
            if g.startswith(prefix):
                r.append(reason)
        if r:
            why[g] = "; ".join(r)
    return why


def oriented(sam, names):
    """sequence of each wanted genome in rCRS orientation, bases other than ACGT written as the rCRS base"""
    out, fixed, removed = {}, 0, 0
    for line in open(sam):
        if line[0] == "@":
            continue
        f = line.split("\t", 11)
        if f[0] not in names or int(f[1]) & (4 | 256 | 2048):
            continue
        seq = list(f[9])
        r, q = int(f[3]) - 1, 0
        if re.search("[^ACGT]", f[9]):
            for n, op in CIGAR.findall(f[5]):
                n = int(n)
                if op == "M":
                    for i in range(n):
                        if seq[q + i] not in "ACGT":
                            seq[q + i] = ref[r + i]
                            fixed += 1
                    r, q = r + n, q + n
                elif op in "IS":
                    for i in range(n):
                        if seq[q + i] not in "ACGT":
                            seq[q + i] = ""
                            removed += 1
                    q += n
                elif op == "D":
                    r += n
        out[f[0]] = "".join(seq)
    return out, fixed, removed


entries, members, seen = [], defaultdict(list), {}


def add(name, seq):
    if seq in seen:
        members[seen[seq]].append(name)
        return False
    seen[seq] = name
    entries.append((name, seq))
    return True


add("chrM", ref)
n_hg, name = 0, None
for line in open(W / "mitoleaf.fa"):
    if line[0] == ">":
        name = line[1:].strip()
    else:
        n_hg += add(name, line.strip())
stats = pd.read_csv(W / "scan/sequence_stats.tsv", sep="\t").set_index("genome")
why = screen(W / "scan", stats)
keep = set(stats.index) - set(why)
seqs, fixed, removed = oriented(W / "scan/aln.sam", keep)
n_gb = sum(add(g, seqs[g]) for g in sorted(seqs))
why_p = screen(W / "scan_people")
rcrs = ref[:3106] + "N" + ref[3106:]
n_pp = 0
for sample, calls in people.read_calls(CALLS).items():
    if sample in why_p:
        continue
    n_pp += add(("HGDP_" if sample.startswith("HGDP") else "1KG_") + sample,
                people.sequence(calls, rcrs, length_indels=False))
with open(OUT, "w") as f:
    for name, seq in entries:
        f.write(f">{name}\n{seq}{seq[:TAIL]}\n")
with open(OUT.parent / "left_out.tsv", "w") as f:
    f.write("sequence\treason\n")
    for g, r in sorted({**why, **why_p}.items()):
        f.write(f"{g}\t{r}\n")
with open(OUT.parent / "members.tsv", "w") as f:
    f.write("entry\tidentical_sequences\n")
    for k, v in members.items():
        f.write(f"{k}\t{' '.join(v)}\n")
print("entries written:", len(entries), "(chrM 1, mitoLEAF haplogroups", n_hg, ", GenBank genomes", n_gb,
      ", re-called people", n_pp, ")")
print("GenBank: on the list", len(stats), ", left out", len(why), ", kept", len(seqs), ", of those distinct and new", n_gb)
print("   bases written as the rCRS base:", fixed, "; removed (in an insertion or clipped):", removed)
print("re-called people: left out", len(why_p), ", distinct and new", n_pp)
print("identical sequences folded into an entry:", sum(len(v) for v in members.values()))
c = Counter()
for r in why.values():
    for part in r.split("; "):
        c[re.sub(r"\d+", "N", part)[:60]] += 1
print("reasons (GenBank):", c.most_common(8))
