#!/usr/bin/env python3
"""One sequence per mitoLEAF haplogroup: rCRS with the haplogroup's substitutions,
written without the N at 3107; identical sequences once, under the first name.
Insertions and deletions in a motif are not applied. Positions mitoLEAF marks as
ambiguous (IUPAC) keep the rCRS base, unless the code excludes it; then the
alphabetically last base of the code is used.

usage: haplogroups.py hgmotifs.json ref.fa > mitoleaf.fa
"""
import json
import re
import sys

IUPAC = {"R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def substitutions(motif, rcrs):
    toks = set()
    for t in motif.split():
        m = re.fullmatch(r"(\d+)([A-Za-z])", t)
        if not m:
            continue
        pos, s = int(m[1]), m[2]
        if s in "ACGT":
            toks.add(t)
        elif s in IUPAC and rcrs[pos - 1] not in IUPAC[s]:
            toks.add(f"{pos}{sorted(IUPAC[s])[-1]}")
    return toks


ref = "".join(l.strip() for l in open(sys.argv[2]) if l[0] != ">")
rcrs = ref[:3106] + "N" + ref[3106:]
seen = set()
for h, motif in json.load(open(sys.argv[1])).items():
    seq = list(rcrs)
    for t in substitutions(motif, rcrs):
        seq[int(t[:-1]) - 1] = t[-1]
    s = "".join(seq[:3106] + seq[3107:])
    if s not in seen:
        seen.add(s)
        print(">mitoLEAF_%s\n%s" % (re.sub(r"[^A-Za-z0-9]", "_", h), s))
