#!/usr/bin/env python3
"""The re-called people as whole sequences, all calls applied, for scan_humans.py.

usage: people_sequences.py people_calls.tsv.gz ref.fa > singles.fa
"""
import sys

import people

ref = "".join(l.strip() for l in open(sys.argv[2]) if l[0] != ">")
rcrs = ref[:3106] + "N" + ref[3106:]
for sample, calls in people.read_calls(sys.argv[1]).items():
    print(f">{sample}\n{people.sequence(calls, rcrs)}")
