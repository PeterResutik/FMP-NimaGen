"""The shared and the separate frame for the C-stretches (the frame spec).

Only 16180-16193 and 300-315 have two frames that differ; 57-60 is laid from
the left like them but has no separate variant. Everywhere else, including
the other complex regions, the general rule in notation.py applies.

Shared frame: a molecule's bases are laid on the region from its first
position onwards, one base per position. A non-C base after the leading run,
where the interrupt would be, anchors the frame: the two sides are laid
separately, with differences left of it at the base before the interrupt and
right of it at the region's last position. Missing bases are deletions at the
3' end of their side, extra bases insertions there. The base closest to the
interrupt position anchors only when that reading needs fewer changes than
without it; each frame decides this on its own reading, and on a tie the
separate frame anchors while the shared frame does not.

Separate frame: a molecule that starts with an unbroken leading run is read
as run lengths instead. A's missing from the leading run are deletions at its
end (A16183-), extra A's an insertion after it, and the C's are laid from the
first C position. A shift at the boundary (A16183C) thus becomes a shorter
leading run plus a longer C-run (A16183- -16193.1C). A leading run broken by
another base (AAACA...) keeps its shared-frame reading.
"""
import re
from collections import Counter, defaultdict
from dataclasses import dataclass

import notation

IUPAC = {frozenset(k): v for k, v in {"AG": "R", "CT": "Y", "AC": "M", "GT": "K", "CG": "S", "AT": "W"}.items()}


@dataclass(frozen=True)
class Region:
    first: int
    last: int
    lead: str = None        # base of the leading run, if any
    lead_end: int = None    # last position of the leading run
    interrupt: int = None   # position of the non-C base inside the C-stretch, if any


REGIONS = {
    "16180-16193": Region(16180, 16193, "A", 16183, 16189),  # A4 C5 T C4
    "300-315": Region(300, 315, "A", 302, 310),              # A3 C7 T C5
    "57-60": Region(57, 60),                                 # T4, laid from the left
}


def _lay(bases, first, n):
    """bases on n positions from first: {pos: base}, deleted positions, bases left over."""
    return {first + i: b for i, b in enumerate(bases[:n])}, [first + i for i in range(len(bases), n)], bases[n:]


def _layout_from_left(region, bases, anchor):
    """Shared frame: laid from the left; with an anchor (index of the interrupt base)
    the two sides are laid separately."""
    if anchor is None:
        placed, deleted, extra = _lay(bases, region.first, region.last - region.first + 1)
        return placed, deleted, {region.last: extra}
    pl, dl, il = _lay(bases[:anchor], region.first, region.interrupt - region.first)
    pr, dr, ir = _lay(bases[anchor + 1:], region.interrupt + 1, region.last - region.interrupt)
    return {**pl, region.interrupt: bases[anchor], **pr}, dl + dr, {region.interrupt - 1: il, region.last: ir}


def _count_changes(layout_, reference):
    placed, deleted, inserted = layout_
    return (sum(b != reference[p - 1] for p, b in placed.items()) + len(deleted)
            + sum(len(seq) for seq in inserted.values()))


def _anchor(region, bases, reference, read):
    """Index of the base that anchors the frame at the interrupt, or None. The candidate
    is the non-C base after the leading run closest to the interrupt in the frame's
    reading (`read`); of two equally close, the one needing fewer changes. It anchors
    when that reading needs fewer changes than the reading without it. On a tie the
    shared frame keeps the reading without, as mitoLEAF writes it (AAAACCCCTCCCCC is
    C16188T T16189C); the separate frame anchors, so the T stays and the run lengths
    change (C16188- -16193.1C)."""
    if not region.interrupt:
        return None
    candidates = [i for i in range(len(_head(region, bases)), len(bases)) if bases[i] != "C"]
    if not candidates:
        return None
    a = _leading_run(region, bases) if read is _separate else None
    if a is None:
        position = lambda i: region.first + i          # laid from the region's first position
    else:
        position = lambda i: region.lead_end + 1 + i - a  # the C's laid from the first C position
    changes = lambda i: _count_changes(read(region, bases, i), reference)
    idx = min(candidates, key=lambda i: (abs(position(i) - region.interrupt), changes(i)))
    anchored, laid = changes(idx), changes(None)
    return idx if anchored < laid or (anchored == laid and a is not None) else None


def _head(region, bases):
    """The leading run, with at most one other base between leading bases (AAACA)."""
    if not region.lead:
        return ""
    m = re.match(rf"{region.lead}+(?:[^{region.lead}]{region.lead}+)?", bases)
    return m.group(0) if m else ""


def _leading_run(region, bases):
    """Length of the unbroken leading run, or None when it is broken (AAACA...)."""
    head = _head(region, bases)
    return None if re.search(rf"[^{region.lead}]", head) else len(head)


def _separate(region, bases, anchor):
    a = _leading_run(region, bases) if region.lead else None
    if a is None or (anchor is not None and anchor < a):
        return _layout_from_left(region, bases, anchor)
    lead_len = region.lead_end - region.first + 1
    placed = {region.first + i: region.lead for i in range(min(a, lead_len))}
    deleted = [region.first + i for i in range(a, lead_len)]
    inserted = {region.lead_end: region.lead * max(0, a - lead_len)}
    if anchor is None:
        pc, dc, extra = _lay(bases[a:], region.lead_end + 1, region.last - region.lead_end)
        placed.update(pc)
        deleted += dc
        inserted[region.last] = extra
    else:
        pl, dl, il = _lay(bases[a:anchor], region.lead_end + 1, region.interrupt - region.lead_end - 1)
        pr, dr, ir = _lay(bases[anchor + 1:], region.interrupt + 1, region.last - region.interrupt)
        placed.update({**pl, region.interrupt: bases[anchor], **pr})
        deleted += dl + dr
        inserted[region.interrupt - 1], inserted[region.last] = il, ir
    return placed, deleted, inserted


def layout(region, bases, frame, reference):
    """(placed {pos: base}, deleted [pos], inserted {anchor: bases}) of one molecule's
    region bases in the 'shared' or 'separate' frame."""
    read = _separate if frame == "separate" else _layout_from_left
    return read(region, bases, _anchor(region, bases, reference, read))


def labels(region, bases, frame, reference):
    """Major labels of one molecule, e.g. ['A16183-', 'T16189C', '-16193.1C']."""
    placed, deleted, inserted = layout(region, bases, frame, reference)
    out = [f"{reference[p - 1]}{p}{b}" for p, b in placed.items() if b != reference[p - 1]]
    out += [f"{reference[p - 1]}{p}-" for p in deleted]
    out += [f"-{a}.{k}{b}" for a, seq in inserted.items() for k, b in enumerate(seq, 1)]
    return sorted(out, key=notation.label_position)


def rows(region, molecules, coverage, frame, reference, min_vf=5.0, lh_thresh=10.0):
    """Report rows [(label, percent)] for the molecules [(region bases, reads)] of one
    sample, out of `coverage` reads. Substitutions use min_vf, the caller's floor
    (--min_vf_FDS for FDSTOOLS molecules; minor as IUPAC), length changes lh_thresh
    (--lh_thresh; minor in lowercase). In the shared frame a C on a leading-run
    position is a boundary shift and uses lh_thresh too, so both frames report the
    same molecules."""
    lh_floor, lh_ceiling = min(lh_thresh, 100 - lh_thresh), max(lh_thresh, 100 - lh_thresh)
    bases_at, deleted_at, inserted_at = defaultdict(Counter), Counter(), defaultdict(Counter)
    for bases, reads in molecules:
        placed, deleted, inserted = layout(region, bases, frame, reference)
        for p, b in placed.items():
            bases_at[p][b] += reads
        for p in deleted:
            deleted_at[p] += reads
        for anchor, seq in inserted.items():
            for k, b in enumerate(seq, 1):
                inserted_at[(anchor, k)][b] += reads
    out = []
    for p in range(region.first, region.last + 1):
        ref = reference[p - 1]
        shares = {b: 100 * n / coverage for b, n in bases_at[p].items()}
        for b, share in sorted(shares.items(), key=lambda x: -x[1]):
            if b == ref:
                continue
            boundary = region.lead and b == "C" and p <= region.lead_end
            floor, ceiling = (lh_floor, lh_ceiling) if boundary else (min_vf, 100 - min_vf)
            if share >= ceiling:
                out.append((f"{ref}{p}{b}", share))
            elif share >= floor:
                out.append((f"{ref}{p}{IUPAC[frozenset((ref, b))]}", share))
        share = 100 * deleted_at[p] / coverage
        if share >= lh_ceiling:
            out.append((f"{ref}{p}-", share))
        elif share >= lh_floor:
            out.append((f"{ref}{p}{ref.lower()}", share))
    for (anchor, k), counts in sorted(inserted_at.items()):
        share = 100 * sum(counts.values()) / coverage
        b = counts.most_common(1)[0][0]
        if share >= lh_ceiling:
            out.append((f"-{anchor}.{k}{b}", share))
        elif share >= lh_floor:
            out.append((f"-{anchor}.{k}{b.lower()}", share))
    return sorted(out, key=lambda r: notation.label_position(r[0]))
