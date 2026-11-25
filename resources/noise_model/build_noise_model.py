#!/usr/bin/env python3
"""A noise model for the depth rule (--noise_model) from control samples: per control, how many wrong bases reach each
share of the reads, every 0.05% from 0.1 to 2% (with at least 2 reads).

usage: build_noise_model.py --library LIBRARY --reference rCRS.fasta [--mask POSITIONS] CONTROL.sast.csv ... > MODEL.tsv

Each CONTROL.sast.csv is one control's FDSTOOLS sequences per amplicon (p12_fdstools/<sample>/<sample>.sast.csv) from a
run of the pipeline with --minimum 1, so that sequences seen once are kept; the control is named by its file name. In
each amplicon the most common sequence is the control's own, and every substitution another sequence has more or
lacks is a wrong base, with the reads of the sequences that carry it. A position's reads are those of every amplicon
that reads it. Left out: insertions and deletions, the C-stretches (57-60, 300-315, 16180-16193) and the control's own
heteroplasmies (--mask; by default those of 2800M: 1282, 2846, 8078, 10000, 12630, 14847, 15759).

A share of 0.1% needs 2 reads in 2,000, so a control with fewer reads per amplicon undercounts the low shares.
The rows of several model files are pooled: a lab can add its controls to the shipped model or use only its own.
"""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import depth_rule  # noqa: E402

MASK_2800M = (1282, 2846, 8078, 10000, 12630, 14847, 15759)
SUBSTITUTION = re.compile(r"^(no )?([ACGT])(\d+)([ACGT])$")
BASES = "ACGT"


def read_reference(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))


def position_of(label):
    return int(re.search(r"\d+", label).group())


def amplicon_sequences(sast_path):
    """{amplicon: [(sequence, reads)]} of a FDSTOOLS samplestats file, most reads first."""
    out, header = {}, None
    for line in open(sast_path):
        fields = line.rstrip("\n").split("\t")
        if header is None:
            header = fields
            continue
        row = dict(zip(header, fields))
        if row["sequence"] in ("Other sequences", "No data"):
            continue
        out.setdefault(row["marker"], []).append((row["sequence"], int(float(row["total"]))))
    for seqs in out.values():
        seqs.sort(key=lambda x: -x[1])
    return out


def count_control(sast_path, segments, reference, masked):
    """(median reads per amplicon, wrong bases reaching each share of depth_rule.SHARES)"""
    reads_at = np.zeros(depth_rule.GENOME + 1)
    wrong = {}                                    # (position, base) -> reads
    own = {}                                      # position -> the control's own base where it is not rCRS
    totals = []
    sequences = amplicon_sequences(sast_path)
    for amplicon, parts in segments.items():
        seqs = sequences.get(amplicon, [])
        total = sum(k for _s, k in seqs)
        totals.append(total)
        for a, b in parts:
            reads_at[a:b + 1] += total
        if not seqs:
            continue
        truth = set(seqs[0][0].split()) - {"REF"}
        for seq, k in seqs:
            variants = set(seq.split()) - {"REF"}
            errors = [x for x in variants - truth] + ["no " + x for x in truth - variants]
            for e in errors:
                m = SUBSTITUTION.match(e)
                if not m or not depth_rule.counted(int(m[3]), masked):
                    continue
                no, ref, p, alt = m.groups()
                p = int(p)
                base = ref if no else alt
                wrong[(p, base)] = wrong.get((p, base), 0) + k
                if no:
                    own[p] = alt
    k_e, n_e = [], []
    for p in range(1, depth_rule.GENOME + 1):
        if reference[p - 1] not in BASES or not depth_rule.counted(p, masked) or reads_at[p] == 0:
            continue
        for base in BASES:
            if base != own.get(p, reference[p - 1]):
                k_e.append(wrong.get((p, base), 0))
                n_e.append(reads_at[p])
    k_e, n_e = np.array(k_e), np.array(n_e)
    counts = [int((k_e >= np.maximum(2, np.ceil(s * n_e - 1e-9))).sum()) for s in depth_rule.SHARES]
    return int(np.median(totals)), counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("sast", nargs="+", help="a control's p12 sast.csv (pipeline run with --minimum 1)")
    parser.add_argument("--library", required=True, help="FDSTOOLS library file (resources/fdstools/)")
    parser.add_argument("--reference", required=True, help="rCRS FASTA (resources/rCRS/)")
    parser.add_argument("--mask", default=",".join(map(str, MASK_2800M)),
                        help="the control's own heteroplasmies, comma-separated (default: 2800M's)")
    args = parser.parse_args()
    masked = sorted({int(p) for p in args.mask.split(",") if p.strip()})
    segments = depth_rule.library_segments(args.library)
    reference = read_reference(args.reference)
    print("# FMP-NimaGen noise model for the depth rule (resources/scripts/depth_rule.py): per control, the wrong "
          "bases with at least 2 reads that reach each share (%)")
    print("# built with build_noise_model.py from " + ", ".join(Path(p).name for p in args.sast))
    print("\t".join(["control", "median_reads_per_amplicon", "masked"] + depth_rule.SHARE_COLUMNS))
    for path in args.sast:
        name = Path(path).name.removesuffix(".sast.csv")
        median, counts = count_control(path, segments, reference, set(masked))
        if median < 2000:
            print(f"{name}: {median} reads per amplicon (median); below 2,000 the shares under "
                  f"{math.ceil(2 / max(median, 1) * 1e4) / 100:g}% are undercounted", file=sys.stderr)
        print("\t".join([name, str(median), ",".join(map(str, masked))] + [str(c) for c in counts]))


if __name__ == "__main__":
    main()
