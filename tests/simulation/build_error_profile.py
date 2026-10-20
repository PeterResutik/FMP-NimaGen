#!/usr/bin/env python3
"""Error profile of real NimaGen reads, for simulate_reads.py --error-profile.

For every source (a sequencing run: the samples given with the same --source) and amplicon it keeps a pool of real
read pairs, each read as its quality string, the substitution errors it carries (cycle, base read) and the
insertion/deletion errors it carries (cycle, "D", number of bases skipped; or cycle, "I", bases inserted).
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
error. An insertion or deletion is the sample's own when at least half of those reads carry it there (placed at the
leftmost position of its homopolymer on rCRS); otherwise it is an error. Not counted as errors (the read keeps its
qualities there): positions where the sample's WGS shows another base, an insertion or a deletion on at least
--mask-af of its reads (an insertion or deletion masks its homopolymer) or --mask-positions lists a known
heteroplasmy (biology, not error), the C-stretch regions (their PCR stutter cannot be told apart from a person's
length heteroplasmy read by read, so it is not simulated), the wobble bases of both of the amplicon's primers (read
at the start and, by R1 or R2 reading through, at the end), and cycles bwa soft-clipped. An insertion or deletion
touching one of these positions or its homopolymer is not counted either.
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


def norm_del(rcrs, p, k):
    """leftmost rCRS position of a deletion of k bases starting at p (shifted along its homopolymer or repeat)"""
    while p > 1 and p + k - 1 <= TRUE_MT_LENGTH and rcrs[p - 2] == rcrs[p + k - 2]:
        p -= 1
    return p


def norm_ins(rcrs, p, s):
    """leftmost rCRS position (and inserted bases) of an insertion s before position p"""
    while p > 1 and rcrs[p - 2] == s[-1]:
        s, p = s[-1] + s[:-1], p - 1
    return p, s


def indel_span(rcrs, kind, p, what):
    """rCRS positions an insertion or deletion at p touches: its bases, one on each side and its homopolymer"""
    k = what if kind == "D" else len(what)
    span = set(range(p - 1, p + k + 1))
    base, q = (rcrs[p - 1] if kind == "D" else what[0]), p
    while q <= TRUE_MT_LENGTH and rcrs[q - 1] == base:
        span.add(q)
        q += 1
    return {(x - 1) % TRUE_MT_LENGTH + 1 for x in span}


def wgs_mask(vcf, min_af, rcrs):
    """rCRS positions where the WGS has a substitution, insertion or deletion allele on at least min_af of the reads
    (mitopy VCF); an insertion or deletion masks the positions it touches."""
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
                if len(ref) != len(alt):
                    k = abs(len(ref) - len(alt))
                    if len(ref) > len(alt):
                        mask |= indel_span(rcrs, "D", norm_del(rcrs, pos + 1, k), k)
                    else:
                        mask |= indel_span(rcrs, "I", *norm_ins(rcrs, pos + 1, alt[1:]))
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
    rcrs = "".join(l.strip() for l in open(args.reference) if not l.startswith(">"))[:TRUE_MT_LENGTH].upper()
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
        mask = (wgs_mask(vcfs[name], args.mask_af, rcrs) if name in vcfs else set()) | extra_mask.get(name, set())
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
            # insertions and deletions: (cycle, "D", bases skipped or "I", bases inserted, rCRS key); a deletion is
            # recorded at the cycle of the base that follows it in the read's own direction
            indels, qp, rp, seq = [], 0, r.reference_start, r.query_sequence
            for op, k in r.cigartuples:
                if op in (0, 7, 8):
                    qp, rp = qp + k, rp + k
                elif op == 2:
                    cyc = (n - qp if r.is_reverse else qp) - off
                    if cyc > 0:
                        p = rp % TRUE_MT_LENGTH + 1
                        indels.append((cyc, "D", k, ("D", norm_del(rcrs, p, k), k)))
                    rp += k
                elif op == 1:
                    cyc = (n - qp - k if r.is_reverse else qp) - off
                    if cyc > 0:
                        ins = seq[qp:qp + k]
                        p, s = norm_ins(rcrs, rp % TRUE_MT_LENGTH + 1, ins)
                        indels.append((cyc, "I", ins.translate(COMP)[::-1] if r.is_reverse else ins, ("I", p, s)))
                    qp += k
                elif op == 4:
                    qp += k
            groups[(amp, d)].append((quals[(r.query_name, d)][off:], cycles, indels))
            pair_reads[(r.query_name, d)] = (amp, len(groups[(amp, d)]) - 1)
            stats[f"primer at offset {off}"] += 1
        for (amp, d), reads in groups.items():
            counts = collections.defaultdict(collections.Counter)
            carriers = collections.Counter()
            for _, cycles, indels in reads:
                for pos, base in cycles.values():
                    counts[pos][base] += 1
                for *_, key in indels:
                    carriers[key] += 1
            major = {k: c.most_common(1)[0][0] for k, c in counts.items()}
            skip = mask | region | wobble[amp]
            # the sample's own insertions and deletions, and those touching masked positions, are not errors
            not_error = set()
            for key in carriers:
                before = (key[1] - 2) % TRUE_MT_LENGTH + 1  # the base before it, read by every read across it
                if carriers[key] >= 0.5 * max(sum(counts[before].values()), carriers[key]) or \
                        indel_span(rcrs, *key) & skip:
                    not_error.add(key)
            for i, (qual, cycles, indels) in enumerate(reads):
                errors = [(cyc, base) for cyc, (pos, base) in sorted(cycles.items())
                          if base != major[pos] and pos not in mask and pos not in region
                          and pos not in wobble[amp]]
                indel_errors = [(cyc, kind, what) for cyc, kind, what, key in indels if key not in not_error]
                reads[i] = (qual, errors, indel_errors)
                stats["reads"] += 1
                stats["errors"] += len(errors)
                stats["indel errors"] += len(indel_errors)
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
        profile["pools"][key] = reads  # [((R1 qualities, R1 errors, R1 indels), (R2 ...))]
    with gzip.open(args.out, "wb") as f:
        pickle.dump(profile, f)
    print(dict(stats))
    print(f"pools: {len(profile['pools'])}; per base: substitution errors {stats['errors'] / max(stats['bases'], 1):.2e}, "
          f"insertion/deletion errors {stats['indel errors'] / max(stats['bases'], 1):.2e}")


if __name__ == "__main__":
    main()
