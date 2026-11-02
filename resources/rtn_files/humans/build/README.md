# Human mitogenomes for rtn, built for NimaGen

rtn compares every read with a set of human mitogenomes and sets a read to
MAPQ 0 (a possible NUMT) when it is far from all of them. A person whose mtDNA
is far from every genome in the set loses real reads, so the set should cover
as many lineages as possible without containing anything NUMT-like.

`../humans_NimaGen_2026-10-30.fa.bz2` (the pipeline's default) is built here from
public sources only. It has 58,346 entries:

| Entries | Source |
|---|---|
| 1 | rCRS (`chrM`) |
| 6,403 | mitoLEAF haplogroups: rCRS with the haplogroup's motif (see below) |
| 49,283 | GenBank genomes on mitoLEAF's list that pass the screen below |
| 2,659 | 1000 Genomes and HGDP people, re-called from their whole-genome data |

A sequence that occurs more than once is written once (17,925 identical
sequences are folded into an entry, listed in `members.tsv`). Every entry is
followed by its own first 100 bases, so reads across the origin still align in
one piece (the NimaGen amplicon over the origin ends at 53).

The earlier list, `../humans_NimaGen_2026-10-19.fa.bz2` (the default up to
pipeline version 0.4.0; uncompressed md5 2c2f66804cec51acf5c138f2ca991ade), is
kept. It was built from the same sources, with haplogroup sequences of the
substitutions only and the rCRS base for two-state marks and missing bases.

## Sources

| Source | Version | Used for |
|---|---|---|
| [mitoLEAF](https://github.com/forensicgenomics/mitoLeaf) `docs/data/hgmotifs.json` | 1.6.1, commit 25b6076, md5 39521b01741850ece285b4806715bbb7 | haplogroup sequences; each genome's haplogroup in the screen |
| mitoLEAF `docs/data/tree.json` | same commit, md5 402b206e62759e2671cdf91ca86af85a | the parent of each haplogroup (for the haplogroup sequences) |
| mitoLEAF `docs/data/profiles.csv` | same commit (tables of 2026-05-26), md5 f7f00e93575fc70e7052330dfc259437 | the 69,714 GenBank accessions (rows with source `NCBI`) |
| GenBank | the accession versions listed in `profiles.csv`, downloaded 2026-10-02 (again 2026-10-20, same file) | genome sequences |
| 1000 Genomes Project, 30x high-coverage CRAMs (Byrska-Bishop et al. 2022) | 2,504 people | `people_calls.tsv.gz` |
| HGDP high-coverage CRAMs (Bergström et al. 2020) | 828 people | `people_calls.tsv.gz` |
| `../../numts/Calabrese_Dayama_Smart_Numts_modified.fa` (from rtn) | as in this repository | the NUMT check |

The people's chrM reads were re-called with mitopy (GATK Mutect2 in
mitochondria mode). A person's sequence is rCRS with their majority calls
(allele fraction 0.5 or more, FILTER `PASS` or `blacklisted_site`);
`people_calls.tsv.gz` holds these calls (`people_calls.py` made it from the
mitopy output). Insertions and deletions in the length regions are not used,
because whole-genome calls are unreliable there, so those stretches have the
rCRS length.

## The screen

A GenBank genome or person is left out (`left_out.tsv`, with the reason) when:

1. it does not align to rCRS in one piece;
2. the insert of a NimaGen amplicon (the stretch between its primers, which
   rtn judges in a primer-trimmed read) is closer to a NUMT than to its own
   mitoLEAF haplogroup by 3 or more differences (`amplicon_check.py`);
3. it has more than 10 unconfirmed substitutions, or 4 or more of them within
   50 bases. A substitution is confirmed when a genome outside the genome's own
   accession series (same prefix, number in the same block of 500) carries it
   too; the genomes counted are those of mitoLEAF's list and the people. A
   variant that two independent genomes share is a lineage variant, not an
   error, even where mitoLEAF's tree does not resolve the lineage;
4. GenBank only: it has more than 5 bases other than A, C, G, T (one to five,
   in practice N, are filled from the sequence of the genome's haplogroup; the
   rCRS base where that haplogroup has a deletion there);
5. GenBank only: it belongs to submission MN628108-MN628167 (Zhang 2019), where
   27 of 69 genomes carry a stretch identical to a NUMT on chromosome 1 or 3.

Left out: 3,185 GenBank genomes (3,028 for ambiguous bases, 76 for the NUMT
check, 69 for the submission, 24 and 12 for rule 3, 8 not in one piece; some
for more than one reason) and no person.

## The haplogroup sequences

Each mitoLEAF haplogroup is written as one sequence: rCRS with the haplogroup's
motif applied (`haplogroups.py`).

- Substitutions (73G), deletions (8281-) and insertions (573.1C) are applied as
  written.
- Lowercase marks are written as present: the optional insertion 315.1c that
  mitoLEAF adds to almost every haplogroup (almost all people carry 315.1C; the
  genome sequences do not record it reliably), and the few others (60.1t,
  745.1t, 960.1c, 573.1c, 16193c and so on).
- A mark that allows two or more states (an IUPAC code, such as 16519Y or the
  16181M 16182M 16183M that mitoLEAF adds to almost every haplogroup) gets the
  state most of the haplogroup's genomes carry: the GenBank genomes and people
  that pass the screen, with the haplogroup `scan_humans.py` gave them. A
  haplogroup without genomes takes the genomes of its nearest ancestral clade
  (the ancestor and all its descendants). A tie goes to the parent's clade,
  counting only the tied states, and further up while it stays tied; rCRS
  only if the whole tree is tied. A genome's ambiguous bases are not counted.

`haplogroup_marks.tsv` lists every such decision: the haplogroup, the mark,
where it was decided (own genomes, an ancestral clade, a tie broken by a
clade), the states of the genomes counted, and the state written.

## Rebuilding

```
bash build.sh workdir
```

`build.sh` needs bwa, samtools and the pipeline's conda environment. It
downloads the three mitoLEAF tables (and checks their md5) and the 69,714 GenBank
genomes (one to two hours; finished batches are kept if it is interrupted), then
runs the steps below. The uncompressed result has md5
9493a1df2a9b5d99481b82ed54ea6bc0; copy `workdir/humans_NimaGen_2026-10-30.fa.bz2`
to `..` and `left_out.tsv`, `members.tsv` and `haplogroup_marks.tsv` here.

| Script | What it does |
|---|---|
| `fetch_genbank.py` | downloads GenBank sequences by accession (NCBI E-utilities) |
| `genbank_singles.py` | one line per genome, and per genome its length and number of ambiguous bases |
| `people_sequences.py`, `people.py` | the people's sequences from `people_calls.tsv.gz` |
| `scan_humans.py` | aligns each genome to rCRS, matches it to the mitoLEAF haplogroup that explains most of its differences, and lists the differences the haplogroup does not explain |
| `amplicon_check.py` | cuts every genome into the inserts of the NimaGen amplicons and flags inserts closer to a NUMT than to the genome's own haplogroup |
| `haplogroups.py` | the rule for the haplogroup sequences (used by `build_humans.py`) |
| `build_humans.py` | applies the screen, builds the haplogroup sequences from the genomes that pass it, and writes the file, `left_out.tsv`, `members.tsv` and `haplogroup_marks.tsv` |

## Licences

- mitoLEAF is licensed under the Mozilla Public License 2.0; the haplogroup
  sequences built from `hgmotifs.json` are covered by that licence.
- GenBank sequences: NCBI places no restrictions on their use; some submitters
  may claim rights in the sequences they submitted.
- 1000 Genomes Project and HGDP data are public without access restrictions;
  please cite Byrska-Bishop et al. 2022 and Bergström et al. 2020.
