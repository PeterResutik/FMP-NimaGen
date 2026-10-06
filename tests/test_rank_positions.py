"""rank_positions.py: positions ranked by branches with a substitution in mitoLEAF's tree."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "resources" / "rtn_files" / "positions"))
import rank_positions as rp  # noqa: E402

TREE = """<?xml version="1.0" encoding="UTF-8"?>
<Node Id="mt-MRCA" HG="73G 263G">
  <Node Id="A" HG="152C 16519C 309.1C 16189Y">
    <Node Id="B" HG="T152 16189C 8281-"/>
  </Node>
  <Node Id="C" HG="152C 195C 73G"/>
</Node>"""


def test_forward_and_back_substitutions_count_and_the_root_does_not():
    c = rp.branch_counts(TREE)
    assert c[152] == 3        # 152C twice, T152 (back) once
    assert c[73] == 1         # only branch C; the root's motif is not a branch
    assert 263 not in c


def test_insertions_deletions_and_ambiguity_codes_do_not_count():
    c = rp.branch_counts(TREE)
    assert c[16189] == 1      # 16189C counts, 16189Y does not
    assert 309 not in c and 8281 not in c


def test_hypervariable_positions_first_then_by_count_and_position():
    order = rp.ranking(rp.branch_counts(TREE))
    assert order[:3] == [16519, 16182, 16183]
    assert order[3:] == [152, 73, 195, 16189]


def test_list_in_the_repository_matches_the_ranking():
    rows = [l.rstrip("\n").split("\t") for l in open(Path(rp.__file__).parent / "recurrent_positions.tsv")][1:]
    assert len(rows) == 50
    assert [int(r[0]) for r in rows] == list(range(1, 51))
    assert [int(r[1]) for r in rows[:5]] == [16519, 16182, 16183, 152, 16311]
