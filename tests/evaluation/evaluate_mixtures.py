#!/usr/bin/env python3
"""Observed against expected for read mixtures of two real samples (mix_reads.py).

usage: evaluate_mixtures.py MIX_OUT PURE_OUT --out PREFIX [--min-vf 2] [--frame separate] [--reference rCRS.fasta]

MIX_OUT: pipeline output of the mixtures (reports and p09 BAMs, --publish_bams true); a mixture is a sample named
mix-<A>-<B>-m<share of B>, optionally followed by _<anything> (mix-HG03295-HG00117-m5_L001).
PURE_OUT: output of the same pipeline version with both samples of every pair on their own (<A>_*_<frame>_frame.xlsx).
--min-vf: the detection threshold of the reports in MIX_OUT, in %; a site whose minor share reaches it is expected
as a heteroplasmy (IUPAC code), otherwise as the major base.

Sites: positions where the two samples' major substitutions differ. Left out: positions
where either sample has a heteroplasmy, a deletion or an insertion of its own, and the
regions 300-315 and 16180-16193 (written per molecule; listed for review in
OUT_PREFIX.length.tsv instead).

Expected: the share of reads from each sample at the position, counted from the read-name
tags in the BAM of p09 (the reads after the NUMT filter); the frequency expected for the
non-rCRS base(s) follows from it. FDSTOOLS counts its reads per amplicon from p08, the same
reads with their primers, so its expected share is about the same.
"""
import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

import openpyxl

REPO = Path(__file__).resolve().parents[2]
TRUE_LENGTH = 16569
MIXTURE = re.compile(r"mix-([^-]+)-([^-]+)-[mB](\d+(?:p\d+)?)(?:_.*)?")
LENGTH_REGIONS = ((300, 315), (16180, 16193))
CODES = {frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("AC"): "M", frozenset("GT"): "K",
         frozenset("CG"): "S", frozenset("AT"): "W", frozenset("CGT"): "B", frozenset("AGT"): "D",
         frozenset("ACT"): "H", frozenset("ACG"): "V"}
SUBSTITUTION = re.compile(r"^([ACGT])(\d+)([ACGTRYMKSWBDHVN])$")
POSITION = re.compile(r"(\d+)")


def report(path):
    rows = [[c.value for c in r] for r in openpyxl.load_workbook(path).active.iter_rows()]
    header = rows[0]
    return [dict(zip(header, r)) for r in rows[1:] if r and r[0] is not None]


def report_of(out, sample, frame):
    found = glob.glob(os.path.join(out, "p13_merged_variants_xlsx", f"{sample}_*{frame}_frame.xlsx"))
    if len(found) != 1:
        sys.exit(f"no single report for {sample} in {out}: {found}")
    return report(found[0])


def position(label):
    m = POSITION.search(str(label))
    return int(m.group(1)) if m else None


def in_length_region(pos):
    return any(lo <= pos <= hi for lo, hi in LENGTH_REGIONS)


def own_calls(rows):
    """(major base per position, positions the sample makes unclear) from a sample's own report."""
    bases, unclear = {}, set()
    for r in rows:
        label = str(r["FMP"])
        if label == "LOW":
            continue
        m = SUBSTITUTION.match(label)
        if m and m.group(3) in "ACGT":
            bases[int(m.group(2))] = m.group(3)
        elif position(label) is not None:
            unclear.add(position(label))
    return bases, unclear


def depth_by_source(bam):
    """{position: (reads of A, reads of B)} from the read-name tags; the appended copy of
    1-53 is added to 1-53."""
    depth = defaultdict(lambda: [0, 0])
    for k, tag in enumerate(("mixA_", "mixB_")):
        view = subprocess.Popen(["samtools", "view", "-u", "-e", f'qname=~"^{tag}"', bam], stdout=subprocess.PIPE)
        out = subprocess.run(["samtools", "depth", "-a", "-Q", "0", "-q", "0", "-"], stdin=view.stdout,
                             capture_output=True, text=True, check=True).stdout
        view.wait()
        for line in out.splitlines():
            _, pos, n = line.split("\t")
            pos = int(pos)
            pos = pos - TRUE_LENGTH if pos > TRUE_LENGTH else pos
            depth[pos][k] += int(n)
    return depth


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mix_out")
    p.add_argument("pure_out")
    p.add_argument("--out", required=True, help="prefix of the four output tables")
    p.add_argument("--min-vf", type=float, default=2, help="detection threshold of the mixture reports, in %%")
    p.add_argument("--frame", choices=["separate", "shared"], default="separate")
    p.add_argument("--reference", default=REPO / "resources/rCRS/rCRS_NimaGen.fasta")
    args = p.parse_args()
    mix_out, pure_out, prefix, min_vf = args.mix_out, args.pure_out, args.out, args.min_vf / 100
    rcrs = "".join(l.strip() for l in open(args.reference) if not l.startswith(">"))[:TRUE_LENGTH]

    site_rows, extra_rows, length_rows, summary = [], [], [], []
    suffix = f"_{args.frame}_frame.xlsx"
    for xlsx in sorted(glob.glob(os.path.join(mix_out, "p13_merged_variants_xlsx", "mix-*" + suffix))):
        sample = os.path.basename(xlsx)[: -len(suffix)]
        m = MIXTURE.fullmatch(sample)
        if not m:
            sys.exit(f"not a mixture name (mix-<A>-<B>-m<share>): {sample}")
        a, b, nominal = m[1], m[2], float(m[3].replace("p", "."))
        nominal = int(nominal) if nominal.is_integer() else nominal
        pure_a, pure_b = own_calls(report_of(pure_out, a, args.frame)), own_calls(report_of(pure_out, b, args.frame))
        mix = report(xlsx)
        by_pos = defaultdict(list)
        for r in mix:
            if str(r["FMP"]) != "LOW" and position(r["FMP"]) is not None:
                by_pos[position(r["FMP"])].append(r)

        bam = glob.glob(os.path.join(mix_out, "p09_filtered_numts_bam_for_mutect2", sample, "*.bam"))[0]
        depth = depth_by_source(bam)

        (bases_a, unclear_a), (bases_b, unclear_b) = pure_a, pure_b
        sites = sorted(p for p in set(bases_a) | set(bases_b)
                       if bases_a.get(p, rcrs[p - 1]) != bases_b.get(p, rcrs[p - 1])
                       and p not in unclear_a | unclear_b and not in_length_region(p))
        errors = []
        counts = defaultdict(int)
        for p in sites:
            ref = rcrs[p - 1]
            base_a, base_b = bases_a.get(p, ref), bases_b.get(p, ref)
            n_a, n_b = depth[p]
            if n_a + n_b == 0:
                counts["no reads"] += 1
                continue
            share_b = n_b / (n_a + n_b)
            minor = min(share_b, 1 - share_b)
            expected_vf = (0 if base_a == ref else 1 - share_b) + (0 if base_b == ref else share_b)
            if minor >= min_vf:
                expected = f"{ref}{p}{CODES[frozenset((base_a, base_b))]}"
            else:
                major = base_a if share_b < 0.5 else base_b
                expected = f"{ref}{p}{major}" if major != ref else None
            rows = [r for r in by_pos.get(p, []) if SUBSTITUTION.match(str(r["FMP"]))]
            observed = rows[0] if rows else {}
            label = observed.get("FMP")
            if expected is None:
                status = "not expected, found" if label else "not expected, absent"
            elif label == expected:
                status = "found"
            elif label is None:
                status = "missed"
            else:
                status = "other label"
            counts[status] += 1
            if minor >= min_vf:
                counts["heteroplasmies expected"] += 1
                counts["heteroplasmies found"] += status == "found"
            vf_fds = number(observed.get("vf_FDS"))
            if status == "found" and minor >= min_vf and vf_fds is not None:
                errors.append(vf_fds - 100 * expected_vf)
            site_rows.append([sample, a, b, nominal, p, ref, base_a, base_b, n_a, n_b, round(100 * share_b, 2),
                              round(100 * expected_vf, 2), expected or "", label or "",
                              vf_fds if vf_fds is not None else "", status])

        explained = set(sites) | unclear_a | unclear_b | set(bases_a) | set(bases_b)
        for r in mix:
            p = position(r["FMP"])
            if str(r["FMP"]) == "LOW" or p is None:
                continue
            if in_length_region(p):
                length_rows.append([sample, nominal, r["FMP"], r.get("vf_FDS") or ""])
            elif p not in explained:
                extra_rows.append([sample, nominal, r["FMP"], r.get("vf_FDS") or "", r.get("marker") or ""])
        n_extra = sum(1 for e in extra_rows if e[0] == sample)
        med = round(median(abs(x) for x in errors), 2) if errors else ""
        bias = round(median(errors), 2) if errors else ""
        summary.append([sample, a, b, nominal, len(sites), counts["heteroplasmies expected"],
                        counts["heteroplasmies found"], counts["found"], counts["missed"], counts["other label"], counts["not expected, found"], med, bias, n_extra])

    def write(name, header, rows):
        with open(f"{prefix}.{name}.tsv", "w") as f:
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")

    write("sites", ["sample", "major", "minor", "minor_pct", "pos", "rCRS", "base_major", "base_minor", "reads_major",
                    "reads_minor", "share_minor_pct", "expected_vf_pct", "expected", "FMP", "vf_FDS", "status"],
          site_rows)
    write("extra", ["sample", "minor_pct", "FMP", "vf_FDS", "marker"], extra_rows)
    write("length", ["sample", "minor_pct", "FMP", "vf_FDS"], length_rows)
    write("summary", ["sample", "major", "minor", "minor_pct", "sites", "het_expected", "het_found", "found", "missed",
                      "other_label", "found_not_expected", "abs_err_FDS", "bias_FDS", "extra_rows"], summary)


if __name__ == "__main__":
    main()
