# Primers and amplicons

The kit (NimaGen IDS-MTC, whole mtDNA) amplifies 100 overlapping amplicons in
two pools (pA, pB). The bases under one amplicon's primers lie in the insert of
a neighbouring amplicon, so every position is read without primer bases
somewhere.

| File | Used by | What it holds |
|---|---|---|
| `left_primers.fasta` | `p04`, `p06` | left oligos, 5' to 3' |
| `right_primers_rc.fasta` | `p04`, `p06` | right oligos as on the forward strand (reverse complement) |
| `../fdstools/mtNG_library_file.txt` | `p12`, `p13` | each amplicon's insert, with its primers as flanks |
| `../amplicon_bed/amplicons_bed.txt` | `p08`, `p09` | one middle position per amplicon, for read depth |
| `../amplicon_bed/NimaGen_primers.bed` | IGV only | primer sites and wobble positions |

Positions are rCRS. `rCRS_NimaGen.fasta` repeats positions 1-53 after 16569, so
that amplicon 100, which crosses the origin, maps in one piece. The BED draws
100's right primer on that copy (16603-16622).

## Primer ends

- 5' ends: every primer begins where the reads begin, in data from both
  chemistries, apart from the exceptions below.
- 3' ends: a primer longer than these files say would leave primer bases in
  the reads, counted as sample bases. To check, we used samples that differ
  from the rCRS close to a primer's end: inside a primer, that primer's reads
  show the primer's base whatever the sample has; past its end, they show the
  sample's base. For the 39 primers where such a sample exists, the reads show
  the sample's base right after the end given here, so no primer is longer
  than in these files. In the last one to three bases of a primer the reads
  also often show the sample's base: the primer may really end there, or the
  polymerase replaced a mismatched primer end with the sample's base. Either
  way, cutting these bases loses nothing, as a neighbouring amplicon reads them.

## Primers at a wrong site

Some primers also bind elsewhere in the mtDNA: reads that start at 3318 begin
with amplicon 083's left primer (home 13644), and the same happens at 12976
(025), 13367 (077) and 16250 (006). bwa clips the part of the primer that does
not fit there and p04 removes it; the part that fits stays in the read, but it
comes from the oligo, not the sample. cutadapt finds no primer at the start of
such a read and drops it, so these reads do not reach the callers.

## Wobbles

A letter other than A, C, G or T (IUPAC code) marks a position where the oligo
pool carries more than one base; the primer part of the reads shows the
mixture in nearly every sample. The primer files have 70 such positions. The
BED marks 64 of them, and the reads show all 64. In the earlier chemistry the
second base is weaker at most of them and absent at some.

The other six (2758, 5773, 8020, 8200, 14560, 15927) lie 2-4 bases from a
primer's 3' end. There the mixture test does not work: a sample is copied
mainly by the oligo that matches it, and the polymerase corrects a mismatch at
the 3' end. Only samples with the other base can tell, and few of ours have
it. With 2758A the amplicon nearly drops out in both chemistries, and with
15927A in the current one, so the kit has no oligo with that base. With 14560A
it drops out in the earlier chemistry but not in the current one, which
probably has the oligo. No sample here carries 5773A, 8020A or 8200C. The codes
stay until NimaGen confirms; a code only lets cutadapt accept either base.

## Exceptions

- **001**: the left oligo begins 5 bases before position 1 (CGATG,
  16565-16569). bwa soft-clips them and `p04` removes them, so the primer file
  and the library flank begin at position 1.
- **045**: the current left oligo covers 7086-7112 (its 3' end is not known
  to the base; the longer of two possible ends is used). The earlier
  chemistry's left oligo covers 7107-7130 (`045_pA_7132-7310_l_short`, primer
  file only). FDSTOOLS reads 045 from 7131, after the earlier oligo, so that
  it reads both chemistries; there 7113-7130 come from amplicon 044.
- **046**: the kit has two left oligos, one from 7273 with C at 7274 and one
  from 7272 with T at 7274. The primer file lists the first, with Y at 7274;
  `p04` cuts the extra first base of the second as read-through. The BED shows
  both.
- **100**: the insert runs from 16471 across the origin to 33.

The coordinates in the primer names are labels; in most amplicons they differ
from the insert by one or two bases.
