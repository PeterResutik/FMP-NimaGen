#!/usr/bin/env python3
"""Symmetric length-heteroplasmy floor/ceiling around a single --lh_thresh
value - a reporting-threshold convention (below floor: not reported,
floor-ceiling: minor/lowercase, at/above ceiling: major/uppercase/dash),
not reference structure. Split out of repeat_regions.py (2026-09-23,
user: "repeat_regions.py contains a lot of functions that are not
related to repeat regions") into its own module, since this genuinely
serves a different purpose than that file's reference-sequence math.

Used by process_fdstools_output_improved_better.py, merge_fdstools_
mutect2_improved.py, call_repeat_regions.py (percentage scale, 0-100),
and process_mutect2_output_improved.py (fraction scale, 0-1, via
lh_bounds - its own apply_snp/apply_insertion/apply_deletion/
finalize_output_table use this convention throughout, args.lh_thresh/100
at the CLI boundary).
"""


def lh_bounds_pct(threshold_pct):
    """Symmetric floor/ceiling (percentage scale, 0-100) around a single
    length-heteroplasmy threshold, e.g. threshold_pct=10 -> (10, 90).
    Accepts either side (10 or 90) and always returns (floor, ceiling)
    with floor <= ceiling."""
    return min(threshold_pct, 100 - threshold_pct), max(threshold_pct, 100 - threshold_pct)


def lh_bounds(threshold):
    """Fraction-scale (0-1) equivalent of lh_bounds_pct above, e.g.
    threshold=0.10 -> (0.10, 0.90). Kept as its own explicitly-scaled
    entry point, a thin wrapper around the one real implementation,
    rather than forcing every fraction-scale call site onto the
    percentage scale for no functional benefit."""
    floor_pct, ceiling_pct = lh_bounds_pct(threshold * 100)
    return floor_pct / 100, ceiling_pct / 100
