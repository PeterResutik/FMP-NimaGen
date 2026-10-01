"""IUPAC codes for the bases present at one position, shared by both callers' scripts,
the frames and the merge: two bases R Y M K S W, three bases V H D B, all four N."""

CODES = {frozenset(bases): code for bases, code in {
    "AG": "R", "CT": "Y", "AC": "M", "GT": "K", "CG": "S", "AT": "W",
    "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B", "ACGT": "N"}.items()}
BASES = {code: set(bases) for bases, code in CODES.items()}
