# Test data

`2800M-DVJ49_S1_L001_R1_001.fastq.gz` and `_R2_` are read pairs of the Promega
2800M positive control of a MiSeq run (flow cell DVJ49, 2 × 226 cycles), for
`-profile test` and the pipeline test (`tests/main.nf.test`): 22,448 of the
control's 142,731 read pairs, about 100 reads per amplicon. Kept are 1 in 6 (by
the MD5 of the read name) of the pairs whose merged read the pipeline maps, and
of the primer dimers. Left out are the pairs that FLASH cannot merge, and the
pairs whose merged read carries, at a position where 2800M has another base,
the base of another sample of the run (1.4% of the merged reads). The reads are
renamed (`2800M-DVJ49.<n>`); their bases and qualities are those of the run.

`humans/humans_test.fa.bz2` holds 227 of the human mitogenomes of
`resources/rtn_files/humans/humans_NimaGen_2026-10-30.fa.bz2`: rCRS, mitoLEAF's
haplogroups with names of one or two characters (the top of every clade), and
H1c with all below it (2800M is H1c). Its index builds in seconds, where the
full file takes about 20 minutes; with it, rtn compares the test reads with
fewer genomes than in a real run.
