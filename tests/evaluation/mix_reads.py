#!/usr/bin/env python3
"""Mix the read pairs of two real samples at an exact share in every amplicon.

usage: mix_reads.py A_R1_001.fastq.gz B_R1_001.fastq.gz share_B_percent out_prefix [seed]

Out: <out_prefix>_R1_001.fastq.gz, <out_prefix>_R2_001.fastq.gz and <out_prefix>_mix.tsv (pairs taken and available
per amplicon). Name the prefix like a sample, e.g. mix-<A>-<B>-m5_L001, for evaluate_mixtures.py.

A pair belongs to the amplicon whose left primer its R1 starts with (the primer's first 15 bases, within the read's
first 6, wobble bases allowed); pairs whose R1 starts with no primer (primer dimers, a read error at the primer
start) form one more group. In every group the most pairs both samples can supply at the share are drawn at random:
round(T * (1 - share)) of A and round(T * share) of B, with T = min(pairs of A / (1 - share), pairs of B / share).
So every amplicon has the same share of B, whatever the two samples' depths per amplicon. Read names get the prefix
mixA_ or mixB_, so the reads of each sample can still be counted after the pipeline. R2 is read alongside R1.
"""
import collections
import gzip
import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "simulation"))
from simulate_reads import IUPAC, REPO, load_primers  # noqa: E402

K = 15


def primer_index(path):
    """{every concrete first-15-bases of a left primer (wobbles expanded): amplicon}"""
    index = {}
    for amp, primers in load_primers(path).items():
        for primer in primers:
            for bases in itertools.product(*(IUPAC.get(c, c) for c in primer[:K])):
                index["".join(bases)] = amp
    return index


def amplicon(seq, index, max_offset=5):
    for off in range(max_offset + 1):
        amp = index.get(seq[off:off + K])
        if amp:
            return amp
    return None


def pairs(r1):
    r2 = r1.replace("_R1_001", "_R2_001")
    with gzip.open(r1, "rt") as f1, gzip.open(r2, "rt") as f2:
        while True:
            a = [f1.readline() for _ in range(4)]
            b = [f2.readline() for _ in range(4)]
            if not a[0]:
                return
            yield a, b


def groups(r1, index):
    """{amplicon or None: [pair numbers]}"""
    out = collections.defaultdict(list)
    for i, (a, _) in enumerate(pairs(r1)):
        out[amplicon(a[1].strip(), index)].append(i)
    return out


def tagged(record, tag):
    return ["@" + tag + record[0][1:]] + record[1:]


def main():
    r1_a, r1_b, pct, prefix = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    share = pct / 100
    index = primer_index(REPO / "resources/primers/left_primers.fasta")
    rng = random.Random(seed)
    group_a, group_b = groups(r1_a, index), groups(r1_b, index)
    keep_a, keep_b, lines = set(), set(), []
    for amp in sorted(set(group_a) | set(group_b), key=lambda a: (a is None, a or "")):
        n_a, n_b = len(group_a.get(amp, [])), len(group_b.get(amp, []))
        total = min(n_a / (1 - share), n_b / share)
        take_a, take_b = min(round(total * (1 - share)), n_a), min(round(total * share), n_b)
        keep_a.update(rng.sample(group_a.get(amp, []), take_a))
        keep_b.update(rng.sample(group_b.get(amp, []), take_b))
        lines.append(f"{amp or 'no_primer'}\t{take_a}\t{n_a}\t{take_b}\t{n_b}")
    with gzip.open(prefix + "_R1_001.fastq.gz", "wt", compresslevel=4) as o1, \
            gzip.open(prefix + "_R2_001.fastq.gz", "wt", compresslevel=4) as o2:
        for r1, keep, tag in ((r1_a, keep_a, "mixA_"), (r1_b, keep_b, "mixB_")):
            for i, (a, b) in enumerate(pairs(r1)):
                if i in keep:
                    o1.writelines(tagged(a, tag))
                    o2.writelines(tagged(b, tag))
    with open(prefix + "_mix.tsv", "w") as f:
        f.write("amplicon\tpairs_A\tof_A\tpairs_B\tof_B\n" + "\n".join(lines) + "\n")
    print(f"{prefix}\tA pairs {len(keep_a)}\tB pairs {len(keep_b)}\tshare of B {pct:g}%")


if __name__ == "__main__":
    main()
