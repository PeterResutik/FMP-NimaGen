#!/usr/bin/env python3
"""Error profile of real NimaGen reads, for simulate_reads.py --error-profile.

For every source (a sequencing run: the samples given with the same --source) and amplicon it keeps a pool of real
read pairs, each read as its quality string and the substitution errors it carries (cycle, base read).
simulate_reads.py gives every simulated haplogroup the noise of one source, and every simulated read pair the
qualities and errors of one real pair of that source and amplicon, so the simulated reads carry the real data's
quality decline along the read, its read-to-read variation, the link between the two reads of a pair and the errors
that recur at fixed cycles in that run. Pools are per run, not per sample: an amplicon with few reads in one sample
would otherwise give the same few real reads to many simulated reads.

Input per sample: the raw FASTQ pair, the pipeline's p02 BAM (raw reads aligned to rCRS) and p09 BAM (reads kept
after the NUMT filter and primer trimming). Only pairs kept in p09 are used, so NUMT reads and reads the pipeline
discards do not count as noise. A read's amplicon is the primer it starts with (first 15 bases, within the read's
first 6; cycles are counted from the primer's first base, as simulate_reads.py writes the reads). The true
base at a position is the majority base of the sample's reads of that amplicon and direction; another base is an
error. Not counted as errors (the read keeps its qualities there): positions where the sample's WGS shows another
base on at least --mask-af of its reads or --mask-positions lists a known heteroplasmy (biology, not error), the
C-stretch regions, the wobble bases of both of the
amplicon's primers (read at the start and, by R1 or R2 reading through, at the end), and
cycles bwa soft-clipped. Insertions and deletions are not recorded.
"""
import argparse, collections, gzip, pickle, random, re
import pysam

IUPAC = {"R": "AG", "Y": "CT", "M": "AC", "K": "GT", "S": "CG", "W": "AT", "B": "CGT", "D": "AGT", "H": "ACT",
         "V": "ACG", "N": "ACGT"}
COMP = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")
TRUE_MT_LENGTH = 16569
REGIONS = [(57, 60), (300, 315), (16180, 16193)]


def load_primers(path):
    """{amplicon id: [primers]}, anchors removed (as simulate_reads.py)."""
    primers, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split("_")[0]
        elif line:
            primers.setdefault(name, []).append(line.strip("^$X").upper())
    return primers


def starts_with(seq, primer, k=15):
    return len(seq) >= k and all(a == b or (b in IUPAC and a in IUPAC[b]) for a, b in zip(seq[:k], primer[:k]))


def find_primer(seq, primers, max_offset=5):
    """(amplicon, offset) of the primer the read starts with, or (None, None)."""
    for off in range(max_offset + 1):
        for amp, ps in primers.items():
            if any(starts_with(seq[off:], s) for s in ps):
                return amp, off
    return None, None


def wobble_positions(reference, left, right):
    """{amplicon: rCRS positions of the wobble bases of its primers}. Left primers end where the insert begins,
    right primers (as on rCRS, right_primers_rc) begin after it; primers located as in simulate_reads.py."""
    rcrs = "".join(l.strip() for l in open(reference) if not l.startswith(">"))[:TRUE_MT_LENGTH].upper()
    circular = rcrs + rcrs[:500]
    regex = lambda p: re.compile("".join(f"[{IUPAC[c]}]" if c in IUPAC else c for c in p))
    out = {}
    for amp in left:
        lm = regex(left[amp][-1]).search(circular)
        rm = regex(right[amp][0]).search(circular, lm.end())
        pos = {(lm.end() - len(p) + i) % TRUE_MT_LENGTH + 1 for p in left[amp] for i, b in enumerate(p) if b in IUPAC}
        pos |= {(rm.start() + i) % TRUE_MT_LENGTH + 1 for p in right[amp] for i, b in enumerate(p) if b in IUPAC}
        out[amp] = pos
    return out


def wgs_mask(vcf, min_af):
    """rCRS positions where the WGS has a substitution allele on at least min_af of the reads (mitopy VCF)."""
    mask = set()
    for line in open(vcf):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        pos, ref, alts = int(c[1]), c[3], c[4].split(",")
        fmt = dict(zip(c[8].split(":"), c[9].split(":")))
        for alt, af in zip(alts, fmt.get("AF", "0").split(",")):
            if float(af) >= min_af:
                for i in range(min(len(ref), len(alt))):
                    if ref[i] != alt[i]:
                        mask.add(pos + i)
    return mask


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", nargs=5, action="append", required=True,
                   metavar=("NAME", "R1", "R2", "P02_BAM", "P09_BAM"))
    p.add_argument("--wgs-vcf", nargs=2, action="append", default=[], metavar=("NAME", "VCF"))
    p.add_argument("--mask-positions", nargs=2, action="append", default=[], metavar=("NAME", "POS,POS"),
                   help="further positions of the sample not counted as errors (known heteroplasmies)")
    p.add_argument("--source", nargs=2, action="append", required=True, metavar=("NAME", "SOURCE"),
                   help="the noise source (sequencing run) of the sample")
    p.add_argument("--mask-af", type=float, default=0.01)
    p.add_argument("--left-primers", required=True)
    p.add_argument("--right-primers", required=True, help="right_primers_rc.fasta")
    p.add_argument("--reference", required=True, help="rCRS_NimaGen.fasta")
    p.add_argument("--pool-size", type=int, default=4000, help="read pairs kept per source and amplicon")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    left, right_rc = load_primers(args.left_primers), load_primers(args.right_primers)
    right = {a: [s.translate(COMP)[::-1] for s in ps] for a, ps in right_rc.items()}
    wobble = wobble_positions(args.reference, left, right_rc)
    region = {x for a, b in REGIONS for x in range(a, b + 1)}
    vcfs = dict(args.wgs_vcf)
    extra_mask = {n: {int(x) for x in v.split(",") if x} for n, v in args.mask_positions}
    source = dict(args.source)

    pools = collections.defaultdict(list)
    stats = collections.Counter()
    for name, fq1, fq2, bam02, bam09 in args.sample:
        kept = {r.query_name for r in pysam.AlignmentFile(bam09) if not r.is_secondary and not r.is_supplementary}
        quals = {}
        for d, fq in (("R1", fq1), ("R2", fq2)):
            with pysam.FastxFile(fq) as f:
                for r in f:
                    if r.name in kept:
                        quals[(r.name, d)] = r.quality
        mask = (wgs_mask(vcfs[name], args.mask_af) if name in vcfs else set()) | extra_mask.get(name, set())
        stats[f"{name} masked positions"] = len(mask)
        groups = collections.defaultdict(list)
        pair_reads = {}
        for r in pysam.AlignmentFile(bam02):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.query_name not in kept:
                continue
            d = "R1" if r.is_read1 else "R2"
            amp, off = find_primer(r.get_forward_sequence(), left if d == "R1" else right)
            if amp is None or (r.query_name, d) not in quals:
                stats["no primer at cycle 1"] += 1
                continue
            n = r.query_length
            # cycle (0-based, in the read's own direction) -> rCRS position (1-based, origin copy folded back)
            cycles = {}
            for qpos, rpos in r.get_aligned_pairs(matches_only=True):
                cyc = (n - 1 - qpos if r.is_reverse else qpos) - off
                if cyc < 0:
                    continue
                base = r.query_sequence[qpos]
                cycles[cyc] = ((rpos % TRUE_MT_LENGTH) + 1, base.translate(COMP) if r.is_reverse else base)
            groups[(amp, d)].append((quals[(r.query_name, d)][off:], cycles))
            pair_reads[(r.query_name, d)] = (amp, len(groups[(amp, d)]) - 1)
            stats[f"primer at offset {off}"] += 1
        for (amp, d), reads in groups.items():
            counts = collections.defaultdict(collections.Counter)
            for _, cycles in reads:
                for pos, base in cycles.values():
                    counts[pos][base] += 1
            major = {k: c.most_common(1)[0][0] for k, c in counts.items()}
            for i, (qual, cycles) in enumerate(reads):
                errors = [(cyc, base) for cyc, (pos, base) in sorted(cycles.items())
                          if base != major[pos] and pos not in mask and pos not in region
                          and pos not in wobble[amp]]
                reads[i] = (qual, errors)
                stats["reads"] += 1
                stats["errors"] += len(errors)
                stats["bases"] += len(qual)
        n_pairs = 0
        for (qname, d), (amp, i) in pair_reads.items():
            if d == "R1" and pair_reads.get((qname, "R2"), (None,))[0] == amp:
                pools[(source[name], amp)].append((groups[(amp, "R1")][i], groups[(amp, "R2")][pair_reads[(qname, "R2")][1]]))
                n_pairs += 1
        stats["pairs"] += n_pairs
        print(f"{name} ({source[name]}): {n_pairs} read pairs", flush=True)

    rng = random.Random(args.seed)
    profile = {"pools": {}, "meta": {"sources": sorted(set(source.values())), "samples": source, "mask_af": args.mask_af,
                                     "regions": REGIONS, "stats": dict(stats)}}
    for key, reads in sorted(pools.items()):
        if len(reads) > args.pool_size:
            reads = rng.sample(reads, args.pool_size)
        profile["pools"][key] = reads  # [((R1 qualities, R1 errors), (R2 qualities, R2 errors))]
    with gzip.open(args.out, "wb") as f:
        pickle.dump(profile, f)
    print(dict(stats))
    print(f"pools: {len(profile['pools'])}; error rate per base {stats['errors'] / max(stats['bases'], 1):.2e}")


if __name__ == "__main__":
    main()
