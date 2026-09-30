#!/usr/bin/env python3
"""Simulate NimaGen amplicon reads for mitoLEAF haplogroups.

Each haplogroup's genome is rCRS with its mitoLEAF motif applied. For every
amplicon, a read pair is sequenced from both ends as in real MiSeq data: R1
reads --read-length bases from the left primer, R2 from the right primer, and
both stop at the amplicon end when the amplicon is shorter. Primer
sequences come from the primer files, as they would after PCR. Sequencing
errors follow per-base qualities; PCR stutter is not simulated.

mitoLEAF marks positions whose state varies within a haplogroup with IUPAC
codes (16519Y: T or C) and lowercase entries (315.1c: C inserted or not;
16193c: C or deleted). With --ambiguous one (default) each such position gets
one of its states, chosen at random but seeded by haplogroup, in all reads.
With --ambiguous het they become 50/50 heteroplasmies on two haplotypes,
reads alternating between them.

Output per haplogroup: <sample>_R1_001.fastq.gz, <sample>_R2_001.fastq.gz
(the pipeline's read pattern) and <sample>_truth.tsv with the variants used.
"""

import argparse
import gzip
import json
import multiprocessing
import random
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TRUE_MT_LENGTH = 16569
IUPAC = {"R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
IUPAC_CODE = {frozenset(bases): code for code, bases in IUPAC.items() if len(bases) == 2}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def load_rcrs(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    return seq[:TRUE_MT_LENGTH].upper()


def load_primers(path):
    """{amplicon id ("001"): [primer sequences]}, anchors (^, $) removed."""
    primers, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split("_")[0]
        elif line:
            primers.setdefault(name, []).append(line.strip("^$").upper())
    return primers


def primer_regex(primer):
    return re.compile("".join(f"[{IUPAC[c]}]" if c in IUPAC else c for c in primer))


def find_inserts(rcrs, left, right):
    """{amplicon id: [rCRS positions between the primers]}, located by
    matching the primers on the circular rCRS (the origin amplicon wraps)."""
    circular = rcrs + rcrs[:500]
    inserts = {}
    for amp in left:
        if amp not in right:
            continue
        lm = primer_regex(left[amp][-1]).search(circular)
        rm = primer_regex(right[amp][0]).search(circular, lm.end()) if lm else None
        if not (lm and rm):
            raise ValueError(f"primers of amplicon {amp} not found on rCRS")
        inserts[amp] = [(i % TRUE_MT_LENGTH) + 1 for i in range(lm.end(), rm.start())]
    return inserts


def apply_motif(rcrs, motif, ambiguous="het", rng=None):
    """Per rCRS position: the base(s) of this haplogroup there ("" if deleted,
    ref + inserted bases for insertions), for two haplotypes. Fixed variants
    are on both. Ambiguous entries (IUPAC, lowercase) are 50/50 with
    ambiguous="het" (variant form on the first haplotype, the other form on
    the second), or one state drawn with rng on both with ambiguous="one".
    Returns (haplotypes, variants used)."""
    # rCRS keeps an N placeholder at 3107 where real genomes have no base
    bases = {p: (b if b in "ACGT" else "") for p, b in enumerate(rcrs, start=1)}
    haplotypes = (bases, dict(bases))
    insertions, used = ({}, {}), []
    for token in motif.split():
        if m := re.fullmatch(r"(\d+)\.(\d+)([ACGTacgt])", token):
            pos, idx, base = int(m[1]), int(m[2]), m[3].upper()
            if m[3].islower() and ambiguous == "one":
                # a second inserted base only follows a first one
                if idx > 1 and idx - 1 not in insertions[0].get(pos, {}):
                    continue
                if rng.random() < 0.5:
                    continue
                token = f"{pos}.{idx}{base}"
            for ins in insertions[:1 if m[3].islower() and ambiguous == "het" else 2]:
                ins.setdefault(pos, {})[idx] = base
            used.append(token)
            continue
        pos = int(re.match(r"\d+", token)[0])
        ref = rcrs[pos - 1]
        if m := re.fullmatch(r"(\d+)([ACGT])", token):
            forms, label = (m[2], m[2]), token
        elif (m := re.fullmatch(r"(\d+)([acgt])", token)) and m[2].upper() == ref:
            # lowercase rCRS base: the base or a deletion
            forms, label = ("", ref), token
            if ambiguous == "one":
                deleted = rng.random() < 0.5
                forms, label = (("", ""), f"{pos}-") if deleted else ((ref, ref), None)
        elif m := re.fullmatch(r"(\d+)([RYMKSWBDHVN])", token):
            alt = next(b for b in IUPAC[m[2]] if b != ref)
            other = ref if ref in IUPAC[m[2]] else next(b for b in IUPAC[m[2]] if b != alt)
            forms, label = (alt, other), f"{pos}{IUPAC_CODE[frozenset((alt, other))]}"
            if ambiguous == "one":
                state = rng.choice(sorted(IUPAC[m[2]]))
                forms, label = (state, state), (f"{pos}{state}" if state != ref else None)
        elif re.fullmatch(r"(\d+)-", token):
            forms, label = ("", ""), token
        else:
            raise ValueError(f"unknown motif notation: {token}")
        for hap, form in zip(haplotypes, forms):
            hap[pos] = form
        if label:
            used.append(label)
    for hap, ins_by_pos in zip(haplotypes, insertions):
        for pos, ins in ins_by_pos.items():
            hap[pos] += "".join(ins[i] for i in sorted(ins))
    return haplotypes, used


def amplicon_template(bases, insert_positions):
    """Insert sequence, including insertions right after the left primer."""
    before = insert_positions[0] - 1 if insert_positions[0] > 1 else TRUE_MT_LENGTH
    return bases[before][1:] + "".join(bases[p] for p in insert_positions)


def resolve_iupac(seq, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in seq)


def sequence_read(molecule, rng, qualities):
    read, qual = [], []
    for base in molecule:
        q = rng.choice(qualities)
        if rng.random() < 10 ** (-q / 10):
            base = rng.choice([b for b in "ACGT" if b != base])
        read.append(base)
        qual.append(chr(q + 33))
    return "".join(read), "".join(qual)


def simulate(haplogroup, motif, sample, args, rcrs, left, right, inserts):
    # Seeded by haplogroup name, so its reads don't depend on batch or order
    rng = random.Random(f"{args.seed}:{haplogroup}")
    # ambiguous states come from their own stream, so read errors don't depend on them
    states = random.Random(f"{args.seed}:{haplogroup}:states")
    haplotypes, used = apply_motif(rcrs, motif, args.ambiguous, states)
    qualities = [33, 34, 35, 36, 36, 37, 37, 37, 38, 38]
    outdir = Path(args.outdir)
    r1 = gzip.open(outdir / f"{sample}_R1_001.fastq.gz", "wt")
    r2 = gzip.open(outdir / f"{sample}_R2_001.fastq.gz", "wt")
    n = 0
    for amp, positions in sorted(inserts.items()):
        templates = [amplicon_template(hap, positions) for hap in haplotypes]
        for i in range(args.reads_per_amplicon):
            # Reads alternate between the haplotypes; primer variants (e.g.
            # amplicon 045) are used in turn within each haplotype
            lp = left[amp][(i // 2) % len(left[amp])]
            molecule = resolve_iupac(lp, rng) + templates[i % 2] + resolve_iupac(right[amp][0], rng)
            n += 1
            name = f"@SIM:1:{sample}:1:{1101 + n // 10000}:{n % 10000}:{amp}"
            s1, q1 = sequence_read(molecule[:args.read_length], rng, qualities)
            s2, q2 = sequence_read(molecule.translate(COMPLEMENT)[::-1][:args.read_length], rng, qualities)
            r1.write(f"{name} 1:N:0:1\n{s1}\n+\n{q1}\n")
            r2.write(f"{name} 2:N:0:1\n{s2}\n+\n{q2}\n")
    r1.close()
    r2.close()
    with open(outdir / f"{sample}_truth.tsv", "w") as truth:
        truth.write(f"haplogroup\t{haplogroup}\nmotif\t{motif}\nambiguous\t{args.ambiguous}\nvariants_used\t{' '.join(used)}\n")
    return n


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--motifs", required=True, help="mitoLEAF hgmotifs.json")
    p.add_argument("--haplogroups", nargs="*", default=[], help="haplogroup names")
    p.add_argument("--haplogroup-list", help="file with one haplogroup name per line")
    p.add_argument("--outdir", required=True)
    p.add_argument("--reads-per-amplicon", type=int, default=100)
    p.add_argument("--read-length", type=int, default=226, help="bases read from each end (current workflow: 226)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--ambiguous", choices=["one", "het"], default="one",
                   help="ambiguous motif entries: one state per haplogroup (one) or 50/50 heteroplasmy (het)")
    p.add_argument("--processes", type=int, default=4)
    p.add_argument("--reference", default=REPO / "resources/rCRS/rCRS_NimaGen.fasta")
    p.add_argument("--left-primers", default=REPO / "resources/primers/left_primers.fasta")
    p.add_argument("--right-primers", default=REPO / "resources/primers/right_primers_rc.fasta")
    args = p.parse_args()

    motifs = json.load(open(args.motifs))
    names = list(args.haplogroups)
    if args.haplogroup_list:
        names += [l.strip() for l in open(args.haplogroup_list) if l.strip()]
    missing = [h for h in names if h not in motifs]
    if missing:
        raise SystemExit(f"not in mitoLEAF: {missing}")

    rcrs = load_rcrs(args.reference)
    left, right = load_primers(args.left_primers), load_primers(args.right_primers)
    inserts = find_inserts(rcrs, left, right)
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    jobs, seen = [], set()
    for hg in names:
        sample = re.sub(r"[^A-Za-z0-9]+", "-", hg).strip("-")
        while sample in seen:
            sample += "x"
        seen.add(sample)
        jobs.append((hg, motifs[hg], sample, args, rcrs, left, right, inserts))
    with multiprocessing.Pool(args.processes) as pool:
        for (hg, _, sample, *_), n in zip(jobs, pool.starmap(simulate, jobs)):
            print(f"{hg} -> {sample}: {n} read pairs over {len(inserts)} amplicons", flush=True)


if __name__ == "__main__":
    main()
