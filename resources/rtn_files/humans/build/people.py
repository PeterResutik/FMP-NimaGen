"""Sequences of the re-called 1000 Genomes and HGDP people from people_calls.tsv.gz
(used by people_sequences.py and build_humans.py)."""
import gzip
from collections import defaultdict

# indels in these regions are length variation; whole-genome calls are unreliable there
LENGTH_REGIONS = [(57, 60), (300, 316), (452, 463), (514, 525), (568, 573), (956, 965),
                  (5895, 5899), (8272, 8289), (16180, 16193), (3105, 3110)]


def read_calls(path):
    """{sample: [(pos, ref, alt, af)]}, samples in the order of the table"""
    calls = defaultdict(list)
    with gzip.open(path, "rt") as f:
        next(f)
        for line in f:
            s, pos, ref, alt, af = line.rstrip("\n").split("\t")
            calls[s].append((int(pos), ref, alt, float(af)))
    return calls


def sequence(calls, rcrs, length_indels=True):
    """rCRS (16,569 bases, N at 3107) with the calls applied from the right, a call that
    overlaps one already applied or does not match rCRS skipped; written without the N.
    length_indels=False leaves out insertions and deletions in LENGTH_REGIONS."""
    seq, left = list(rcrs), len(rcrs) + 1
    for pos, ref, alt, af in sorted(calls, key=lambda c: (-c[0], -c[3])):
        if not length_indels and len(ref) != len(alt) and any(a - 1 <= pos <= b for a, b in LENGTH_REGIONS):
            continue
        if pos + len(ref) - 1 >= left or "".join(rcrs[pos - 1:pos - 1 + len(ref)]) != ref:
            continue
        seq[pos - 1:pos - 1 + len(ref)] = [alt] + [""] * (len(ref) - 1)
        left = pos
    seq[3106] = "" if seq[3106] == "N" else seq[3106]
    return "".join(seq)
