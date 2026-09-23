#!/usr/bin/env python3
"""Derive calls inside reference repeat regions directly from reads.

Inside homopolymer/repeat runs, both callers' own bookkeeping becomes
unreliable (fragmented multiallelic splits whose frequencies don't sum
sensibly, local-reassembly artifacts inventing indel lengths present in
no read, del+sub merges folding distinct events together). Outside such
runs both callers agree with the reads to within ~1-2%, so this module
deliberately only takes over inside reference-derived runs and leaves
everything else untouched.

Calls are derived from LOCAL HAPLOTYPES, not per-position base counts.
Each read's actual observed sequence across the run is extracted and
aligned to the reference run, and the resulting edits are tallied.

That distinction is not cosmetic. An earlier version counted bases per
reference position, restricted to reads with net==0 length change over
the run, on the assumption that "no net length change" meant "read
through cleanly, so its base calls are trustworthy". That assumption is
false when an insertion and a deletion cancel out: HG01799 carries
16182.1C (+1) together with T16189DEL (-1), so its main molecule is
net==0 and bwa-mem aligns it as a gapless 183M - no indel anywhere in
the CIGAR to detect. The bases simply sit one position right of where
they belong, and per-position counting reported the register shift as
three substitutions (A16183C, C16184A, T16189C) that do not exist. All
of that is forensically wrong: EMPOP nomenclature needs the insertion
and the deletion, and the substitution spelling matches no database
entry.

Aligning the whole observed run sequence instead recovers the true
edits, because the shift becomes visible as a gap:

    ref   AAA-ACCCCCTCCCC
    read  AAACACCCCC-CCCC   -> 16182.1C + T16189DEL

Gaps are scored cheaper than mismatches here, which is correct inside a
repeat run (indels are the common event) and is what lets the aligner
prefer one insertion plus one deletion over three substitutions.

Two independent axes, deliberately not collapsed into one number:

  substitutions   - over reads that actually have a base at that
                    position (reads deleted there carry no base call and
                    are excluded from the denominator, not counted as
                    "not the substitution")
  length variants - over all reads spanning the run

Collapsing them is what produces uninterpretable figures: at a position
where most reads carry a deletion (e.g. rCRS 310, the D310 tract), a
single blended percentage answers neither question.

Multi-base indels are reported cumulatively ("at least N") at the run's
rightmost position, matching the forensic 3'-shift convention and the
way FDSTOOLS counts the same events.
"""

import argparse
import re
from collections import Counter, defaultdict

from Bio import Align

from read_evidence import _iter_reads, CIGAR_RE
from reference_utils import load_reference, find_homopolymer_runs, merge_runs
from heteroplasmy_thresholds import lh_bounds_pct

IUPAC_CODES = {
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["A", "C"]): "M",
    frozenset(["G", "T"]): "K",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
}

# Scoring is calibrated against two real, opposite cases, because the
# balance between "spell this as indels" and "spell this as substitutions"
# is exactly what decides whether the output is forensically right:
#
#   HG01799      AAAACCCCCTCCCC -> AAACACCCCCCCCC
#                2 gap events beat 3 mismatches  -> 16182.1C + T16189DEL
#   s25-04643-E2 AAAACCCCCTCCCC -> AACCCCCCCCCCCC
#                3 mismatches beat 3 gap events  -> A16182C A16183C T16189C
#
# Both are confirmed against FDSTOOLS' haplotypes. open_gap=-5 is the only
# band that satisfies both (at -2/-3/-4 the second case collapses into a
# pile of indels). Don't retune these without re-checking both.
_ALIGNER = Align.PairwiseAligner(
    mode="global",
    match_score=2,
    mismatch_score=-3,
    open_gap_score=-5,
    extend_gap_score=-0.5,
    target_end_gap_score=-10,
    query_end_gap_score=-10,
)


def local_haplotypes(bam, chrom, start, end):
    """Observed sequence of every read spanning [start, end], tallied.

    Inserted bases are included in the returned sequence, so a molecule
    carrying an insertion is longer than the reference run - that length
    difference is exactly the signal the alignment step needs. Only reads
    covering both ends are counted, so partial coverage can't masquerade
    as a length change.
    """
    haplotypes = Counter()
    for pos, cigar, seq in _iter_reads(bam, f"{chrom}:{start}-{end}"):
        ref_pos, idx = pos, 0
        chunk = []
        covered_start = covered_end = False
        for length, op in CIGAR_RE.findall(cigar):
            length = int(length)
            if op in ("M", "=", "X"):
                for k in range(length):
                    p = ref_pos + k
                    if p == start:
                        covered_start = True
                    if p == end:
                        covered_end = True
                    if start <= p <= end:
                        chunk.append(seq[idx + k])
                ref_pos += length
                idx += length
            elif op == "I":
                # An insertion sits between reference positions; keep it
                # if the base to its left is inside the run.
                if start <= ref_pos - 1 <= end:
                    chunk.append(seq[idx:idx + length])
                idx += length
            elif op in ("D", "N"):
                # A deletion still *covers* these reference positions - the
                # read spans them, it just has no base there. Missing this
                # silently drops every read whose deletion happens to sit on
                # the run's first or last position, i.e. exactly the reads
                # a deletion call depends on.
                for k in range(length):
                    p = ref_pos + k
                    if p == start:
                        covered_start = True
                    if p == end:
                        covered_end = True
                ref_pos += length
            elif op == "S":
                idx += length
        if covered_start and covered_end:
            haplotypes["".join(chunk).upper()] += 1
    return haplotypes


def decompose(ref_run, observed, start):
    """Align one observed run sequence to a reference run.

    Returns (indels, deleted_positions, bases):
      indels           - ("INS", anchor_pos, inserted_seq) inserted after
                         anchor_pos, and ("DEL", first_pos, deleted_seq)
      deleted_positions- reference positions the read has no base at, so
                         substitution denominators can exclude them
      bases            - {reference position: base the read shows there},
                         which the caller compares against whichever
                         reference it wants to report differences from
    """
    indels = []
    deleted_positions = set()
    bases = {}

    if observed == ref_run:
        for k, base in enumerate(ref_run):
            bases[start + k] = base
        return indels, deleted_positions, bases

    alignment = _ALIGNER.align(ref_run, observed)[0]
    target, query = alignment[0], alignment[1]

    ref_pos = start
    i = 0
    while i < len(target):
        if target[i] == "-":
            run = ""
            while i < len(target) and target[i] == "-":
                run += query[i]
                i += 1
            indels.append(("INS", ref_pos - 1, run))
        elif query[i] == "-":
            run = ""
            first = ref_pos
            while i < len(target) and query[i] == "-":
                run += target[i]
                deleted_positions.add(ref_pos)
                ref_pos += 1
                i += 1
            indels.append(("DEL", first, run))
        else:
            bases[ref_pos] = query[i]
            ref_pos += 1
            i += 1
    return indels, deleted_positions, bases


def structural_decompose(observed, reference, start, end):
    """Decompose one observed run-sequence against the reference's OWN
    leading-run structure, instead of aligning it against a frame
    borrowed from some other molecule (decompose()/_ALIGNER above).

    Why this exists: decompose() needs a reference to align against, and
    the only reference available inside a run is one baked from some
    molecule's own substitutions (typically the sample's dominant one).
    That is provably wrong whenever the sample has more than one
    subpopulation differing on a run-boundary substitution - in
    s26-02989 ~76% of molecules carry A16182C and ~24% don't, and
    EVERY choice of frame is wrong for one group: bake A16182C and the
    24% get a fabricated A-insertion (one gap beats substitute-back-
    plus-gap against that frame); don't bake it and the 76% get a
    fabricated A-deletion (same reasoning, opposite direction). The
    alignment step itself is not malfunctioning in either case - "insert
    an A" genuinely is the minimal edit relative to a 2-A frame when the
    molecule has 3 A's. The frame is what's wrong, not the algorithm.

    This reads each molecule's own structure directly instead of
    borrowing anyone else's frame, so there is nothing to get wrong per
    molecule:

      1. Count the molecule's own leading run (e.g. leading A's, 0-4).
         That count IS the number of reference A-positions still
         unconverted for this molecule - no frame needed to know it.
      2. Whatever's left ("the rest") is read as one block spanning
         start+leading_run through end. Its length compared to the
         reference's own length there gives any real indel, which lands
         exactly on the run's rightmost position by construction (it's
         measured from the end), matching the forensic 3'-shift
         convention without a separate shifting pass.

    Only implemented for the "one leading run, then the remainder of the
    reference span" shape - which is exactly what reference_utils.py's
    merge_runs produces (a run merged with whatever follows across
    single-base gaps), i.e. every region this module is ever given.
    Genuinely wider structures (a second independently-ambiguous run
    later in the span) aren't specifically resolved by this, though
    they're read no worse than plain per-character comparison would.
    A leading-run count that would exceed the reference run's own length
    (an insertion WITHIN the leading run) isn't handled either - treated
    as implausible per se (why would a run insert into itself rather
    than being explained as length change in the neighbouring run,
    which is always the cheaper explanation) and unobserved in every
    validated sample.
    """
    first_base = reference[start - 1].upper()
    first_len = 1
    while start - 1 + first_len < len(reference) and reference[start - 1 + first_len].upper() == first_base:
        first_len += 1
    first_len = min(first_len, end - start + 1)

    matched = 0
    while matched < first_len and matched < len(observed) and observed[matched] == first_base:
        matched += 1
    rest_observed = observed[matched:]
    rest_start = start + first_len
    rest_ref_len = end - rest_start + 1
    rest_ref = "".join(reference[rest_start - 1:end]).upper()

    variants = []
    bases = {}
    deleted_positions = set()

    for i in range(matched):
        bases[start + i] = first_base
    # Unconverted positions get the base "the rest" actually starts
    # with - i.e. whatever this molecule really shows there - reported
    # at the RIGHTMOST unmatched leading-run positions (matches the
    # existing greedy-from-the-left / substitution-on-the-right
    # convention, e.g. 2-of-4 A's converted reports as A16182C+A16183C).
    absorbed_base = rest_observed[0] if rest_observed else first_base
    for i in range(first_len - matched):
        bases[start + matched + i] = absorbed_base

    net = len(rest_observed) - rest_ref_len
    if net == 0:
        for i, ch in enumerate(rest_observed):
            bases[rest_start + i] = ch
    elif net > 0:
        for i in range(rest_ref_len):
            bases[rest_start + i] = rest_observed[i]
        variants.append(("INS", end, rest_observed[rest_ref_len:]))
    else:
        k = -net
        # Everything up to the deletion reads at face value; the
        # rightmost k reference positions are the ones convention places
        # the deletion at, and carry no base call - not "whatever
        # substitution pattern the rest of the region shows", since
        # which physical base is actually missing is exactly what's
        # ambiguous in a homopolymer (matches the module's existing rule
        # that deletion-affected positions are excluded from the
        # substitution axis rather than guessed at).
        for i, ch in enumerate(rest_observed):
            bases[rest_start + i] = ch
        deleted_positions.update(range(end - k + 1, end + 1))
        variants.append(("DEL", end - k + 1, rest_ref[-k:] if k <= len(rest_ref) else rest_ref))

    return variants, deleted_positions, bases


def shift_insertion_right(reference, anchor, segment):
    """Move an insertion to its 3'-most equivalent placement.

    Deliberately NOT the same function as reference_utils.py's own
    shift_insertion_right/shift_deletion_right (used by process_mutect2_
    output_improved.py and process_fdstools_output_improved_better.py) -
    those are VCF-REF/ALT-oriented (a two-phase jump-then-rotate search);
    this one is a simpler single-pass rotation over an already-observed
    read/local-haplotype segment, this file's own different context. Kept
    separate rather than unified since the two algorithms' equivalence
    across multi-base segments was never verified, and this file's own
    output (still a real, standalone tool - see call_boundary_run/
    call_region) wasn't worth that risk for a naming tidy-up alone
    (2026-09-23, user: "pull shared primitives").
    """
    seg = segment
    while anchor < len(reference) and reference[anchor] == seg[0]:
        seg = seg[1:] + seg[0]
        anchor += 1
    return anchor, seg


def shift_deletion_right(reference, first, segment):
    """Move a deletion to its 3'-most equivalent placement."""
    seg = segment
    while first + len(seg) - 1 < len(reference) and reference[first + len(seg) - 1] == seg[0]:
        seg = seg[1:] + seg[0]
        first += 1
    return first, seg


# Regions where process_mutect2_output_improved.py's disable_homopolymer_
# length_calls / merge_fdstools_mutect2_improved.py's mark_disabled_
# mutect2_calls drop/flag Mutect2's own indel/length-axis rows, matching
# the two regions process_fdstools_output_improved_better.py's own
# target_regions covers (2026-09-23, user: "I think we should extend it
# to 303-315 for fdstools and mutect2"). Named for what it now IS - a
# curated, hand-validated allowlist, the Mutect2-side twin of FDSTOOLS'
# own target_regions - not "BOUNDARY_RUN_REGIONS" (renamed 2026-09-23,
# user: "easy to mix up when reading cold"): that name collided with,
# and was easily confused with, find_boundary_run_regions() above - a
# completely different thing, a DYNAMIC, genome-wide detector of every
# boundary-run-SHAPED region the reference happens to contain (9 of
# them), not this file's own curated 2-region validated subset. Each
# entry's "extension" span is where any length-axis row gets dropped;
# "leading" is only set for a TRUE boundary-run region - a leading
# homopolymer run of one base immediately followed by a run of a
# DIFFERENT base (e.g. 16180-16183's A-run into 16184-16193's C-run) -
# where that leading run's own per-position reference-presence also
# needs dropping. chrM 303-315 (the "310" region: 303-309 C's, 310 T
# interrupt, 311-315 C's) has no such leading run - only 3 A's at
# 300-302, below the length-4 homopolymer floor - so its "leading" is
# None; it's still a real homopolymer region Mutect2 struggles with,
# just not boundary-run-shaped (2026-09-23, user: "it's okay that the
# homopolymer stretch is different and doesn't have a leading
# homopolymer stretch").
MUTECT2_TARGET_REGIONS = [
    {"leading": (16180, 16183), "extension": (16184, 16193)},
    {"leading": None, "extension": (303, 315)},
]


def call_boundary_run(bam, chrom, reference, leading_start, leading_end,
                       extension_start, extension_end, min_vf_pct=5.0, lh_thresh_pct=10.0):
    """Call a leading homopolymer run and its immediately-adjacent
    extension run as two independent, position-relative counts, reading
    each read's own local haplotype directly against the true reference -
    no frame is ever baked or borrowed from any other read.

    Named for the region SHAPE it targets (matching find_boundary_run_
    regions above), not the homopolymer_reporting MODE - unlike its own
    FDSTOOLS-side twin, report_separate_frame in process_fdstools_output_
    improved_better.py. This call_boundary_run itself is no longer
    reachable from the automated pipeline at all, only from this file's
    own standalone CLI, since bam_override - the mode that used to call
    it - was removed; kept named after the region shape rather than
    renamed to match report_separate_frame, since it isn't part of that
    mode dispatch (2026-09-23, "pull shared primitives" discussion). See
    report_separate_frame's own docstring for the full design rationale -
    this mirrors it exactly, adapted to raw local-haplotype sequences from
    local_haplotypes() instead of FDSTOOLS' own pre-parsed haplotype
    tokens. Both of call_region's documented known limits are structural
    consequences of needing ONE borrowed alignment frame to spell every
    read; this function never builds one:

    - The LEADING run (e.g. 16180-16183) is reported as independent
      per-position reference-presence counts. Which SPECIFIC position(s)
      lost reference is inherently ambiguous inside a homogeneous run (a
      read missing 2 of 4 A's reads identically regardless of which 2) -
      matching this module's own existing 3'-shift convention (ambiguous
      loss anchored at the RIGHTMOST position(s)), a read's own count of
      LEADING matching characters directly gives the leftmost N positions
      as retained and the rightmost as lost, with no frame needed.

    - The EXTENSION run (e.g. 16184-16193) is reported as one shared
      cumulative "at least k extra bases" count over whatever remains
      after the leading run's own matched prefix - folding in boundary
      insertions (an insertion of the extension base sitting INSIDE the
      leading run, e.g. HG01799's 16182.1C, ends up counted here exactly
      like one anchored at the run's own far end, since it's read as part
      of "the rest" the moment the leading match stops) exactly the same
      way a genuine insertion past the run's own far end is.

    Why this avoids call_region's documented known limits instead of
    repeating structural_decompose's failure: structural_decompose tried
    to build ONE combined per-read description (a bases dict covering
    every position, needed to report substitutions, indels and everything
    else together), and for a read whose leading match falls short of the
    full run, it had to GUESS which base the unmatched leading positions
    "really" showed (its own "absorbed_base" mechanism) - that guess is
    exactly what fabricated HG01799's variants. This function never
    builds that combined description: the leading run and extension run
    are two completely separate tallies, and a position with ambiguous
    reference-loss is simply reported as "reference absent here", with no
    guess about what replaced it - there is nothing left to get wrong.
    """
    haplotypes = local_haplotypes(bam, chrom, leading_start, extension_end)
    total = sum(haplotypes.values())
    if not total:
        return []

    first_base = reference[leading_start - 1].upper()
    first_len = leading_end - leading_start + 1
    extension_base = reference[extension_start - 1].upper()
    extension_ref_len = extension_end - extension_start + 1
    lh_floor, lh_ceiling = lh_bounds_pct(lh_thresh_pct)

    calls = []

    # ---- Leading run: independent per-position reference-presence ----
    # matched is this read's own count of leading_base characters, NOT a
    # consecutive-from-the-left walk: a walk stops the instant it hits an
    # insertion sitting INSIDE the leading run (e.g. HG01799's real BAM
    # sequence is "AAACA...", 3 A's - inserted C - 1 more A - it never
    # reaches 4 consecutive A's even though all 4 reference A's are
    # genuinely present), which is exactly the failure structural_
    # decompose's docstring documents. A read has no explicit tokens
    # telling us where an insertion is anchored the way FDSTOOLS' own
    # haplotype strings do (see report_separate_frame's twin of this
    # function) - simply counting how many leading_base characters exist
    # in the read AT ALL, capped at the reference's own count, sidesteps
    # needing to know where. (This assumes nothing else in the extension
    # run coincidentally introduces an extra leading_base character of
    # its own, e.g. a C->A substitution somewhere in the C-run - not
    # observed in any validated sample, and implausible for this specific
    # locus, but a real, undetected scope limit if it ever occurred.)
    # Everything else in the read - the true extension-run content AND
    # any insertion of the extension base wherever it's anchored - is
    # simply "the read's own length minus however many leading_base
    # characters it has", with no positional bookkeeping needed for that
    # either.
    retained_at = defaultdict(int)
    net_counts = defaultdict(int)
    for observed, n in haplotypes.items():
        matched = min(observed.count(first_base), first_len)
        for i in range(matched):
            retained_at[i] += n
        net_counts[len(observed) - matched - extension_ref_len] += n

    for i in range(first_len):
        pos = leading_start + i
        ref_total = retained_at[i]
        loss_pct = 100.0 * (total - ref_total) / total
        if loss_pct < lh_floor:
            continue
        is_major = loss_pct >= lh_ceiling
        calls.append({
            "label": f"{first_base}{pos}" + ("-" if is_major else first_base.lower()),
            "type": "DEL" if is_major else "LHP",
            "frequency": round(loss_pct, 1),
            "count": total - ref_total,
            "spanning_reads": total,
            "position": pos,
        })

    # ---- Extension run: unified left+right cumulative "at least k" ----
    max_net = max(net_counts) if net_counts else 0
    k = 1
    while k <= max_net:
        reads = sum(n for net, n in net_counts.items() if net >= k)
        if not reads:
            break
        pct = 100.0 * reads / total
        if pct < lh_floor:
            break
        is_major = pct >= lh_ceiling
        base = extension_base if is_major else extension_base.lower()
        calls.append({
            "label": f"-{extension_end}.{k}{base}",
            "type": "INS" if is_major else "LHP",
            "frequency": round(pct, 1), "count": reads,
            "spanning_reads": total, "position": extension_end,
        })
        k += 1

    return sorted(calls, key=lambda c: (c["position"], c["label"]))


def call_region(bam, chrom, reference, start, end, min_vf_pct=5.0, lh_thresh_pct=10.0):
    """Variant calls across one reference repeat run, from local haplotypes."""
    haplotypes = local_haplotypes(bam, chrom, start, end)
    total = sum(haplotypes.values())
    if not total:
        return []

    ref_run = "".join(reference[start - 1:end]).upper()

    # Two references are needed here, built differently, because they
    # answer different questions. Using one for both breaks one of the two
    # real cases this is calibrated against.
    #
    # ALIGNMENT reference - decides how each haplotype is spelled. Bake in
    # the substitutions carried by the single most common haplotype, i.e.
    # the actual dominant molecule, so every other haplotype is read as a
    # difference from it rather than being explained in whatever way scores
    # best on its own (which in s25-04643-E2 put a deletion in the A-run
    # that no other haplotype supports). Its indels are deliberately NOT
    # baked - they are events to report, not backbone.
    #
    # Take the substitutions whether or not the dominant molecule also
    # carries an indel: in HG02389 it has A16183C, T16189C *and* an
    # insertion, and skipping the substitutions there left the insertion
    # unable to shift past 16183, stranding it at 16182.1 instead of
    # 16193.1. In HG01799 the dominant molecule happens to have no
    # substitutions at all, so this bakes nothing and the raw reference is
    # used - which is what that sample needs, since baking T16189C would
    # remove the T and let the aligner spell the register shift as
    # substitutions again, resurrecting the exact bug this module exists to
    # avoid. Never derive this from per-position majority: at HG01799's
    # 16183/16184 the register-shifted bases ARE the majority.
    dominant = max(haplotypes.items(), key=lambda kv: kv[1])[0]
    _dom_indels, _dom_deleted, dom_bases = decompose(ref_run, dominant, start)
    alignment_reference = list(reference)
    for pos, base in dom_bases.items():
        alignment_reference[pos - 1] = base

    # KNOWN LIMIT - a single reference frame cannot spell a run where two
    # subpopulations differ at the same position. In s26-02989 roughly 76%
    # of molecules carry A16182C and 24% do not, and every frame is wrong
    # for one group: bake it and the molecules lacking it get an invented
    # A-insertion at 16181.1 (one gap beats substitute-back-plus-gap);
    # un-bake it and the molecules carrying it get an invented A-deletion
    # (again one gap beats substitution-plus-gap). Both spellings also
    # displace a real 16193.1C insertion. Un-baking was tried and traded
    # one error set for another, so the simpler frame is kept.
    #
    # A structural, count-the-A-run alternative (structural_decompose,
    # further down in this file) was also tried, on the theory that
    # reading each molecule against the reference's OWN structure - rather
    # than a frame borrowed from the dominant molecule - would sidestep
    # this entirely. It does fix s26-02989, but it assumes an indel always
    # sits at or after the leading run's own boundary, which is false for
    # HG01799: that sample's dominant molecule carries 16182.1C - an
    # insertion INSIDE the nominal A-run boundary, between 16182 and the
    # original 16183 - so naive leading-run counting misreads the shifted
    # true 16183 as if it were part of "the rest", and fabricates
    # variants cascading from there. That's not a narrower miss than this
    # module's own known limit; it regressed nearly every sample. Kept
    # here as a documented dead end, not deleted, so it isn't tried again
    # without remembering why it failed.
    ref_run_baked = "".join(alignment_reference[start - 1:end])

    # Targeted fix for the specific failure above, rather than a wholesale
    # replacement: if a molecule's OWN alignment against the dominant frame
    # places an indel INSIDE the leading run (e.g. "insert an A" at
    # 16181.1, when the reference already has 4 A's there), that is a
    # signal the dominant frame's OWN baked substitutions in that run are
    # wrong for this particular molecule - not evidence of a real indel.
    # For s26-02989's 18 minority reads, the dominant frame has both
    # A16182C and A16183C baked; a molecule with only one of those two
    # substitutions gets explained as "the A-run is one base longer" (one
    # gap beats substitute-back-plus-gap against that frame) instead of
    # correctly as "this molecule has A16183C but not A16182C, plus a
    # genuine 16193.1C insertion" - which is what it actually is once
    # re-read against a frame that isn't assuming the disputed position.
    #
    # This is deliberately narrower than structural_decompose: it only
    # discards the LEADING RUN's own baked substitutions, only for
    # molecules whose first attempt already implicated that run, and only
    # after trying the normal (fully baked) frame first - so a real indel
    # that legitimately belongs inside the leading run, like HG01799's
    # 16182.1C, is unaffected (nothing is baked in that sample's A-run to
    # begin with, so this correction is a no-op there).
    first_base = reference[start - 1].upper()
    first_len = 1
    while start - 1 + first_len < len(reference) and reference[start - 1 + first_len].upper() == first_base:
        first_len += 1
    first_len = min(first_len, end - start + 1)
    leading_run = range(start, start + first_len)

    ref_run_no_leading = list(alignment_reference)
    for pos in leading_run:
        ref_run_no_leading[pos - 1] = reference[pos - 1]
    ref_run_no_leading = "".join(ref_run_no_leading[start - 1:end])

    def decompose_molecule(observed):
        indels, deleted, bases = decompose(ref_run_baked, observed, start)
        # The trigger is the indel's own SEGMENT starting with the leading
        # run's base (e.g. an inserted/deleted "A"), not just its anchor
        # POSITION falling in the leading run's numeric range - the two
        # look similar but aren't: a genuine, unrelated C-run event can
        # perfectly well get its naive (pre-shift) anchor placed at e.g.
        # 16182 too, simply because decompose() is free to place a gap
        # anywhere within a long homogeneous C-block, and 16182 reads as
        # 'C' in the baked frame regardless of what it is in raw
        # reference. Checking anchor-in-range alone reroutes those
        # unrelated C-run events too (confirmed: it did, turning a correct
        # C16193c deletion into a fabricated A16181 one) - checking the
        # segment's own base content is what actually isolates the
        # A-run-substitution-ambiguity cases from everything else.
        if any(anchor in leading_run and seg[0] == first_base for _kind, anchor, seg in indels):
            indels, deleted, bases = decompose(ref_run_no_leading, observed, start)
        return indels, deleted, bases

    # What resolves the REMAINING disagreement (once this correction is
    # applied) is knowing which variants the rest of the sample shares,
    # i.e. joint haplotype naming - what FDSTOOLS does against a curated
    # allele set, and what pairwise alignment against any single frame
    # structurally cannot. Treat FDSTOOLS as authoritative over this
    # module when they still disagree in this specific region.

    # SHIFTING reference - decides where an indel is anchored, and must
    # reflect near-universal substitutions even when they are not baked for
    # alignment. In HG01799 the deleted T at 16189 sits between two C-runs
    # that only merge once T16189C (100% of base-called reads) is applied;
    # without this it cannot move and lands at 16189 instead of the 3'-most
    # position 16193 that EMPOP and FDSTOOLS use.
    major_cutoff = 100 - min_vf_pct
    called_here = defaultdict(int)
    base_votes = defaultdict(lambda: defaultdict(int))
    for observed, n in haplotypes.items():
        _i, deleted, bases = decompose_molecule(observed)
        for pos in range(start, end + 1):
            if pos not in deleted:
                called_here[pos] += n
        for pos, base in bases.items():
            base_votes[pos][base] += n
    reference_baked = list(alignment_reference)
    for pos, counts in base_votes.items():
        if not called_here[pos]:
            continue
        base, n = max(counts.items(), key=lambda kv: kv[1])
        if 100.0 * n / called_here[pos] >= major_cutoff:
            reference_baked[pos - 1] = base

    # Pass 2: re-read every haplotype against the backbone. Indels now sit
    # where they belong and shift correctly; substitutions are reported by
    # comparing the bases actually observed against the ORIGINAL reference,
    # so baking changes how things are spelled without hiding any variant.
    lh_floor, lh_ceiling = lh_bounds_pct(lh_thresh_pct)
    sub_counts = defaultdict(int)        # (pos, ref, alt) -> reads
    base_called = defaultdict(int)       # pos -> reads with a base there
    ins_counts = defaultdict(int)        # (anchor, depth, base) -> reads at least this long
    del_counts = defaultdict(int)        # (pos, base) -> reads

    for observed, n in haplotypes.items():
        indels, deleted, bases = decompose_molecule(observed)

        # KNOWN LIMIT, not fixed - substitution-axis exclusion uses
        # decompose()'s NAIVE (pre-shift) deleted-position set, matching
        # FDSTOOLS' own convention: it excludes T16189DEL reads from
        # T16189C's denominator specifically (277/277=100%, HG01799),
        # using each event's own naive anchor, not the 3'-shifted
        # position the LENGTH axis reports it at.
        #
        # The problem: our naive anchor and FDSTOOLS' don't always agree.
        # decompose() is free to place a single-base deletion ANYWHERE
        # within a uniform run with identical alignment score, and in
        # s26-02989 it happened to place several molecules' deletion at
        # naive position 16182 rather than 16189, excluding 45 reads from
        # A16182's denominator when only 8 (whose deletion genuinely
        # reaches back that far) should be excluded. FDSTOOLS apparently
        # anchors each event more specifically (likely via its own
        # independent per-position analysis, not visible to us) and gets
        # this right.
        #
        # Tried switching to the SHIFTED position instead: that fixes the
        # s26-02989 asymmetry but breaks HG01799 the same way in the
        # other direction - T16189DEL's deletion shifts to 16193, so
        # nothing excludes it from 16189 anymore, and (once positions with
        # no naive base are correctly filled from the baked reference
        # rather than silently dropped) T16189C's denominator goes from
        # 277/277 to includes-everything, contradicting FDSTOOLS' own
        # figure. Both attempts are real regressions on a previously-exact
        # case, not narrower misses, so neither is applied. Which naive
        # position "owns" a given deletion isn't recoverable from the
        # observed sequence alone in a homopolymer - it needs the kind of
        # per-event bookkeeping FDSTOOLS has and this module doesn't.
        for pos in range(start, end + 1):
            if pos not in deleted:
                base_called[pos] += n
        for pos, base in bases.items():
            original = reference[pos - 1].upper()
            if base != original:
                sub_counts[(pos, original, base)] += n
        for kind, anchor, seg in indels:
            if kind == "INS":
                anchor, seg = shift_insertion_right(reference_baked, anchor, seg)
                for depth in range(1, len(seg) + 1):
                    ins_counts[(anchor, depth, seg[depth - 1])] += n
            else:
                # Re-read the deleted bases from the shifting reference. The
                # alignment says "a T was deleted at 16189", but once
                # T16189C is baked that position is a C, and deleting a C is
                # what actually lets it shift through the merged run to
                # 16193. Shifting the literal aligned segment instead leaves
                # it stranded at its naive position.
                seg = "".join(reference_baked[anchor - 1:anchor - 1 + len(seg)])
                anchor, seg = shift_deletion_right(reference_baked, anchor, seg)
                for k, base in enumerate(seg):
                    del_counts[(anchor + k, base)] += n

    calls = []

    for (pos, ref_base, alt), n in sub_counts.items():
        denom = base_called[pos]
        if not denom:
            continue
        pct = 100.0 * n / denom
        if pct < min_vf_pct:
            continue
        if pct >= major_cutoff:
            label, vtype = f"{ref_base}{pos}{alt}", "SNP"
        else:
            code = IUPAC_CODES.get(frozenset([ref_base, alt]), f"{ref_base}/{alt}")
            label, vtype = f"{ref_base}{pos}{code}", "PHP"
        calls.append({
            "label": label, "type": vtype, "frequency": round(pct, 1),
            "count": n, "base_called_reads": denom, "position": pos,
        })

    for (anchor, depth, base), n in ins_counts.items():
        pct = 100.0 * n / total
        if pct < lh_floor:
            continue
        is_major = pct >= lh_ceiling
        calls.append({
            "label": f"-{anchor}.{depth}{base if is_major else base.lower()}",
            "type": "INS" if is_major else "LHP",
            "frequency": round(pct, 1), "count": n,
            "spanning_reads": total, "position": anchor,
        })

    for (pos, base), n in del_counts.items():
        pct = 100.0 * n / total
        if pct < lh_floor:
            continue
        is_major = pct >= lh_ceiling
        calls.append({
            "label": f"{base}{pos}" + ("-" if is_major else base.lower()),
            "type": "DEL" if is_major else "LHP",
            "frequency": round(pct, 1), "count": n,
            "spanning_reads": total, "position": pos,
        })

    return sorted(calls, key=lambda c: (c["position"], c["label"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bam")
    parser.add_argument("reference")
    parser.add_argument("--chrom", default="chrM")
    parser.add_argument("--region", help="Single region START-END; default is every reference run")
    parser.add_argument("--min-run", type=int, default=4)
    parser.add_argument("--max-gap", type=int, default=1)
    parser.add_argument("--min_vf", type=float, default=5.0)
    parser.add_argument("--lh_thresh", type=float, default=10.0)
    parser.add_argument("--show-haplotypes", action="store_true",
                        help="Also print the observed local haplotypes behind each region's calls")
    args = parser.parse_args()

    seq = load_reference(args.reference).upper()

    if args.region:
        start, end = (int(x) for x in args.region.split("-"))
        regions = [(start, end)]
    else:
        regions = [(s, e) for s, e, _ in merge_runs(find_homopolymer_runs(seq, args.min_run), args.max_gap)]

    for start, end in regions:
        calls = call_region(args.bam, args.chrom, seq, start, end, args.min_vf, args.lh_thresh)
        if not calls:
            continue
        print(f"=== {args.chrom}:{start}-{end}  ref={seq[start-1:end]} ===")
        if args.show_haplotypes:
            haps = local_haplotypes(args.bam, args.chrom, start, end)
            tot = sum(haps.values())
            for h, c in haps.most_common(8):
                print(f"      {h:<20} {c:>5}  ({100.0*c/tot:.1f}%)")
        for c in calls:
            detail = (f"{c['count']}/{c['base_called_reads']} base-called"
                      if "base_called_reads" in c
                      else f"{c['count']}/{c['spanning_reads']} spanning")
            print(f"  {c['label']:<14} {c['type']:<4} {c['frequency']:>6.1f}%  ({detail})")


if __name__ == "__main__":
    main()
