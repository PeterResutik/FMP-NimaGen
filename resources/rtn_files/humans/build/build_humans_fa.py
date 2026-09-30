#!/usr/bin/env python3
"""Build rtn's humans.fa: the current genomes minus an exclusion list, plus
mitoLEAF haplogroup sequences, each genome written once followed by its own
first --tail bases (so reads across the genome's end still align in one piece).

A haplogroup's sequence is rCRS with its mitoLEAF motif applied and no base for
the N at 3107. Positions mitoLEAF marks as ambiguous (IUPAC, lowercase) keep the
rCRS state; where an IUPAC code excludes the rCRS base, the base other
haplogroups carry more often at that position is used.

usage: build_humans_fa.py humans.fa hgmotifs.json rCRS.fasta exclude.txt add.txt out.fa [--tail 100]
  exclude.txt: genome names to leave out (first column; # comments)
  add.txt:     mitoLEAF haplogroups to add (first column)
Writes out.fa and out.fa.mitoleaf.tsv (label -> haplogroup).
"""
import argparse
import json
import re
from collections import Counter, defaultdict

IUPAC = {"R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def names(path):
    return [l.split()[0] for l in open(path) if l.strip() and not l.startswith("#")]


def read_fasta(path):
    name, parts = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name:
        yield name, "".join(parts)


def haplogroup_sequence(rcrs, motif, fixed_counts):
    bases = {p: b for p, b in enumerate(rcrs, start=1)}
    bases[3107] = ""  # rCRS's N placeholder: real genomes have no base here
    inserted = defaultdict(dict)
    for token in motif.split():
        m = re.fullmatch(r"(\d+)(?:\.(\d+))?([A-Za-z-])", token)
        pos, idx, s = int(m[1]), m[2], m[3]
        if s in IUPAC and rcrs[pos - 1] not in IUPAC[s]:
            s = max(sorted(IUPAC[s]), key=lambda b: fixed_counts[f"{pos}{b}"])
        elif s not in "ACGT-":
            continue  # ambiguous: rCRS state
        if idx:
            inserted[pos][int(idx)] = s
        else:
            bases[pos] = "" if s == "-" else s
    return "".join(bases[p] + "".join(inserted[p][i] for i in sorted(inserted[p]))
                   for p in range(1, len(rcrs) + 1))


def write(out, name, seq, tail):
    seq = seq + seq[:tail]
    out.write(f">{name}\n")
    for i in range(0, len(seq), 60):
        out.write(seq[i:i + 60] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("humans", "motifs", "rcrs", "exclude", "add", "out"):
        ap.add_argument(a)
    ap.add_argument("--tail", type=int, default=100)
    args = ap.parse_args()

    rcrs = "".join(s for _, s in read_fasta(args.rcrs))[:16569]
    motifs = json.load(open(args.motifs))
    fixed_counts = Counter(t for m in motifs.values() for t in m.split() if re.fullmatch(r"\d+[ACGT]", t))
    exclude, add = set(names(args.exclude)), names(args.add)
    missing = [h for h in add if h not in motifs]
    if missing:
        raise SystemExit(f"not in {args.motifs}: {', '.join(missing)}")

    kept, left_out = 0, set()
    with open(args.out, "w") as out:
        for name, seq in read_fasta(args.humans):
            if name in exclude:
                left_out.add(name)
                continue
            half = len(seq) // 2
            if len(seq) % 2 == 0 and seq[:half] == seq[half:]:
                seq = seq[:half]  # written twice in the source
            write(out, name, seq, args.tail)
            kept += 1
        with open(args.out + ".mitoleaf.tsv", "w") as labels:
            labels.write("label\thaplogroup\n")
            used = set()
            for hg in add:
                label = "mitoLEAF_" + re.sub(r"[^A-Za-z0-9]", "_", hg)
                while label in used:
                    label += "_"
                used.add(label)
                write(out, label, haplogroup_sequence(rcrs, motifs[hg], fixed_counts), args.tail)
                labels.write(f"{label}\t{hg}\n")
    if exclude - left_out:
        raise SystemExit(f"excluded genomes not found: {', '.join(sorted(exclude - left_out))}")
    print(f"{kept} genomes kept, {len(left_out)} left out, {len(add)} mitoLEAF haplogroups added")


if __name__ == "__main__":
    main()
