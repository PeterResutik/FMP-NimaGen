# Human mitogenomes for rtn, built for NimaGen

rtn compares every read with a set of human mitogenomes and sets a read to
MAPQ 0 (a possible NUMT) when it is far from all of them. The scripts here
build that set for NimaGen, `humans_NimaGen.fa`, from rtn's own `humans.fa`
(`../humans.fa.bz2`: 43,438 genomes, each written twice):

1. **Four genomes are left out** (`exclude.txt`, with the reason for each): two
   that contain NUMT sequence and would let NUMT reads through, and two that
   are not modern human lineages.
2. **580 mitoLEAF haplogroups are added** (`mitoleaf_added.txt`): those for
   which some NimaGen amplicon has no genome in `humans.fa` matching the
   haplogroup exactly, so that reads of rare lineages are not taken for NUMTs.
   A haplogroup's sequence is rCRS with its motif applied and no base for the N
   at 3107; positions mitoLEAF marks as ambiguous (IUPAC, lowercase) keep the
   rCRS state.
3. **Each genome is written once, followed by its own first 100 bases**, so
   reads across the origin still align in one piece (the NimaGen amplicon over
   the origin ends at 53). This halves the file and the bwa index.

mitoLEAF: <https://github.com/forensicgenomics/mitoLeaf>, version 1.6.1
(commit 25b6076), `docs/data/hgmotifs.json` (md5
39521b01741850ece285b4806715bbb7). mitoLEAF is licensed under the Mozilla
Public License 2.0; the sequences built from it are covered by that licence.

## Rebuilding

```
bunzip2 -k ../humans.fa.bz2
python build_humans_fa.py ../humans.fa hgmotifs.json ../../../rCRS/rCRS_NimaGen.fasta \
    exclude.txt mitoleaf_added.txt humans_NimaGen.fa --tail 100
```

## How the two lists were made

All steps work in one output directory (`out` below) and need bwa, samtools
and the pipeline's conda environment.

| Script | What it does |
|---|---|
| `scan_humans.py humans.fa hgmotifs.json out` | aligns each genome (one copy) to rCRS, writes its differences in mitoLEAF notation, matches it to the haplogroup that explains most of them, and lists the differences the haplogroup does not explain (`scan.tsv`) |
| `numt_check.py out numts.fa` | for genomes with 4 or more such differences within 50 bp, or 25 or more in all, compares that stretch with rCRS and with rtn's NUMT database |
| `tile_check.py out numts.fa` | cuts every genome into 200-bp tiles and flags tiles closer to a NUMT than to the genome's own haplogroup |
| `needed.py out hgmotifs.json mtNG_library_file.txt exclude.txt` | per haplogroup and NimaGen amplicon, differences to the nearest complete genome (`needed.tsv`); `mitoleaf_added.txt` is every haplogroup whose worst amplicon is 1 or more differences away |

`numts.fa` is `../../numts/Calabrese_Dayama_Smart_Numts_modified.fa` and
`mtNG_library_file.txt` is `../../../fdstools/mtNG_library_file.txt`.
