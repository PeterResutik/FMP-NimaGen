"""iupac.py: the codes for the bases present at one position."""
import itertools

import iupac


def test_every_set_of_two_or_more_bases_has_one_code():
    sets = [frozenset(c) for n in (2, 3, 4) for c in itertools.combinations("ACGT", n)]
    assert set(iupac.CODES) == set(sets)
    assert len(set(iupac.CODES.values())) == len(sets) == 11


def test_bases_give_the_code_back():
    for code, bases in iupac.BASES.items():
        assert iupac.CODES[frozenset(bases)] == code
