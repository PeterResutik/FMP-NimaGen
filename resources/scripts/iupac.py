"""IUPAC codes for the bases present at one position, shared by both callers' scripts,
the frames and the merge: two bases R Y M K S W, three bases V H D B, all four N.

A substitution label carries the base itself when it is the only one present (A73G),
otherwise the code of all bases present, rCRS included when present (T16189Y for
T and C, T16189H for T, C and A, T16189M for C and A without T). Its frequency and
reads are given per base other than rCRS: a plain number for one base, "C 40, A 10"
(largest first) for several."""
import re

CODES = {frozenset(bases): code for bases, code in {
    "AG": "R", "CT": "Y", "AC": "M", "GT": "K", "CG": "S", "AT": "W",
    "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B", "ACGT": "N"}.items()}
BASES = {code: set(bases) for bases, code in CODES.items()}

_VALUE = re.compile(r"([ACGT]) ([0-9.]+(?:e-?\d+)?)")


def code(present):
    """What a label carries for the bases present: the base when it is the only one,
    else their code."""
    present = frozenset(present)
    return next(iter(present)) if len(present) == 1 else CODES[present]


def call(ref, shares, min_vf):
    """(code, {base: share}) at one position from the share of every base there
    ({base: percent}, rCRS included): the bases with at least min_vf are present,
    and the label carries the base when it is the only one, else their code. The
    shares of the bases other than rCRS come largest first. None when rCRS is the
    only base present."""
    present = {b for b, s in shares.items() if s >= min_vf}
    others = sorted(present - {ref}, key=lambda b: (-shares[b], b))
    if not others:
        return None
    return code(present), {b: shares[b] for b in others}


def alts(label_code, ref):
    """The bases other than rCRS that a label's base or code stands for."""
    return BASES.get(label_code, {label_code}) - {ref}


def format_values(values, digits=2):
    """A call's vf or reads ({base: value}, bases other than rCRS): a plain number
    for one base, "C 40, A 10" (largest first) for several."""
    if len(values) == 1:
        return round(next(iter(values.values())), digits)
    ordered = sorted(values.items(), key=lambda x: (-x[1], x[0]))
    return ", ".join(f"{b} {round(v, digits):g}" for b, v in ordered)


def parse_values(value, bases):
    """{base: value} for the bases other than rCRS, from format_values' output;
    None when it does not name exactly these bases."""
    if len(bases) == 1:
        try:
            return {next(iter(bases)): float(value)}
        except (TypeError, ValueError):
            pass
    found = {b: float(v) for b, v in _VALUE.findall(str(value))}
    return found if set(found) == set(bases) else None


def scale_values(value, factor, digits=2):
    """format_values' output with every value multiplied by factor (Mutect2's
    fractions to percent); a missing value stays missing."""
    if value is None or (isinstance(value, float) and value != value):
        return value
    try:
        return round(float(value) * factor, digits)
    except (TypeError, ValueError):
        return ", ".join(f"{b} {round(float(v) * factor, digits):g}" for b, v in _VALUE.findall(str(value)))
