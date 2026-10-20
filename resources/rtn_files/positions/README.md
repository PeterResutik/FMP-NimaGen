# Positions rtn should not count

rtn sets a read to MAPQ 0 (a possible NUMT) when it differs from every human
mitogenome in the humans file by about three or more full-quality mismatches. A
person who is not in that file differs from their nearest relative in it mostly
at positions that mutate again and again across the human tree. NUMTs differ
from human mtDNA at positions spread over the whole molecule.

`recurrent_positions.tsv` lists the 50 positions that change most often. In the
copy of the reads that rtn evaluates, these positions are meant to get base
quality 2 (`#`, the lowest value Illumina assigns), so that a difference there
counts about a tenth of a full mismatch. The original qualities stay in the
reads the callers see. This needs an rtn that uses base qualities: rtn as
released ignores them; the fork at <https://github.com/PeterResutik/RtN>
(branch `fixes`) uses them.

## Ranking

First 16519, 16182 and 16183, which are among the fastest-changing substitution
positions in human mtDNA and which PhyloTree and mitoLEAF therefore leave out of
the tree. Then every other position by the number of branches of mitoLEAF's tree
(version 1.6.1, 6,409 haplogroups) with a substitution there, forward or back;
insertions, deletions and ambiguity codes are not counted, ties go by position.
mitoLEAF leaves five disease-associated variants out of its tree (11778 in three
branches, 14484, 14674); this does not affect the ranking.

## Number of positions

Chosen by a rule fixed before any result was seen, on a leave-out of the humans
file (`../humans/humans_NimaGen_2026-10.fa`): a tenth of its genomes and people is
removed at a time (10 folds, 52,843 in all), and one error-free read per removed
entry and amplicon, of the insert only, goes through rtn (with base qualities,
default `-L`). The cost is measured on 663 reads made from rtn's NUMT file and on
the NimaGen amplicons found in GRCh38 and T2T-CHM13 (nuclear copies).

> The smallest N, in steps of 5, where the next 5 positions rescue fewer than 2%
> of the genomes still losing an amplicon, provided that NUMT reads and
> nuclear-copy reads passing rise by no more than 1% over N = 0.

That is N = 50 (`curve.tsv`: up to 45 every step of 5 rescues 2.2% or more; the
step from 50 to 55 rescues none).

| | N = 0 | N = 50 |
|---|---|---|
| leave-out: genomes losing an amplicon, of 52,843 | 158 | 73 |
| NUMT reads passing, of 663 | 319 | 322 |
| nuclear-copy reads passing | 801 | 801 |
| independent test: Mitotree haplogroups losing an amplicon, of 52,953 | 28 | 5 |

The humans file was rebuilt with its NUMT check per amplicon
(`../humans/humans_NimaGen_2026-10-19.fa`). The same rule on the rebuilt file
gives N = 50 again: genomes losing an amplicon 158 at N = 0, 80 at N = 45, 74 at
N = 50 and 55; NUMT reads passing 318 and 321; nuclear-copy reads 801; Mitotree
haplogroups losing an amplicon at N = 50: the same 5.

The Mitotree (FamilyTreeDNA, CC BY-NC-ND 4.0) was used only as a test set, after
N was fixed; nothing from it is in this repository. On 17 real samples the 50
positions change no call.

## Files

| File | Content |
|---|---|
| `recurrent_positions.tsv` | the 50 positions: rank, position, number of branches |
| `rank_positions.py` | the ranking from mitoLEAF's tree |
| `curve.tsv` | the leave-out, N = 0 to 100 in steps of 5 |

## Rebuilding

```
curl -L https://raw.githubusercontent.com/forensicgenomics/mitoLeaf/25b6076/docs/data/mitoLEAF_phm.xml -o mitoLEAF_phm.xml
# md5 b2b53b6c7c5b3754d813d1126a5a5f6a
python rank_positions.py mitoLEAF_phm.xml 50 > recurrent_positions.tsv
```

mitoLEAF is licensed under the Mozilla Public License 2.0.
