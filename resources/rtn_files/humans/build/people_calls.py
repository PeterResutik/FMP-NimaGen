#!/usr/bin/env python3
"""Make people_calls.tsv from the re-called 1000 Genomes and HGDP people: one row per
majority call that a person's sequence uses (mitopy postprocessed VCF, allele
fraction 0.5 or more, FILTER PASS or blacklisted_site, alternative allele of A, C,
G, T only). Only needed to remake the table; build.sh starts from the table.

usage: people_calls.py mitopy_dir > people_calls.tsv
"""
import glob
import os
import re
import sys

MIN_AF = 0.5
OK_FILTERS = {"PASS", "blacklisted_site"}

print("sample\tpos\tref\talt\taf")
for vcf in sorted(glob.glob(os.path.join(sys.argv[1], "samples/*/*_postprocessed.vcf"))):
    sample = os.path.basename(os.path.dirname(vcf))
    for l in open(vcf):
        if l[0] == "#":
            continue
        f = l.rstrip("\n").split("\t")
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        try:
            af = float(fmt["AF"])
        except ValueError:
            continue
        if af >= MIN_AF and set(f[6].split(";")) <= OK_FILTERS and re.fullmatch("[ACGT]+", f[4]):
            print(f"{sample}\t{f[1]}\t{f[3]}\t{f[4]}\t{fmt['AF']}")
