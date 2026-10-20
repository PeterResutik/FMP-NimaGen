# Changelog

All notable changes to this project will be documented in this file.

This project follows [Semantic Versioning](https://semver.org) and uses the [Keep a Changelog](https://keepachangelog.com) format.

---

## [Unreleased]

### Changed
- FLASH counts the mismatches over the whole overlap of a read pair
  (`--max_overlap` 226, the read length; was 140). FLASH also merges overlaps
  longer than `--max_overlap`, but counted their mismatches over the first
  `--max_overlap` bases only; with 2x226 reads nearly every overlap was longer
  than 140, and FLASH warned about it in every run. No overlap is longer than a
  read, so with the read length the whole overlap counts; with longer reads,
  set `--max_overlap` to their length. A few more read pairs are merged, among
  them pairs whose mismatches cluster in one part of the overlap, as in
  amplicon 023; its recurring errors reported by FDSTOOLS alone rise by a few
  tenths of a percent.

## [0.2.0] – 2026-10-16

### Fixed
- The parameter log lists the trimming parameters, `--reads` and `--outdir`
  under the names they are set with (`--quality_cutoff`,
  `--primer_min_overlap`, `--minimum_length`, `--maximum_length`); it showed
  cutadapt's own option names, which the pipeline does not take.
- The parameter log lists `--reference` and `--fdstools_library` under these
  names. It showed them as `mtDNA reference genome` and `--marker_map`, a
  name the pipeline does not take. `--fdstools_library` is now listed under
  FDSTOOLS, which reads it first.
- An interrupted first run no longer leaves an incomplete humans index that
  later runs take as complete. `p01b` checked for the first file bwa index
  writes (`.amb`); it now checks for the last (`.sa`), and unpacks the archive
  under a temporary name, so a cut-off unpacking is redone as well.
- rtn and Mutect2 tasks reserve memory (`--rtn_memory`, `--mutect2_memory`,
  4 GB each), so Nextflow starts no more of them than the computer's memory
  holds. On a laptop one rtn task per core started at once, each holding the
  humans index (about 2 GB), and the run ran out of memory. Mutect2's Java
  memory follows `--mutect2_memory` (it was fixed at 16 GB); its calls are the
  same.
- FDSTOOLS tasks reserve as many cores as they use threads (`--num_threads`, 6,
  at most as many as the computer has), so Nextflow does not start more of
  them than the computer has cores.
  Mutect2's pair HMM runs on 1 thread (`--native_pair_hmm_threads`, was 4):
  more threads saved almost no time on mtDNA, and the calls are the same.
- FDSTOOLS' depth and frequency at a position covered by two amplicons count
  the reads of both. They counted only the amplicons that list the label, so a
  label one amplicon lists and the other does not was measured against one
  amplicon alone: a sequencing error on 2.8% of the reads of amplicon 053 at
  8541 is under 1% of all reads there, as Mutect2 counts it. `marker` and
  `marker_range` now name every amplicon counted, `num_markers` gives their
  number.
- R1 keeps the anchor sequence before the poly-C tract (chrM:278-302) and is
  cut after it (`p02`); before, the anchor was cut off too. 278-302 was then
  read by R2 alone, right after it had read the poly-C tract, where its
  quality drops: a low-quality error at 297 (C instead of A, in up to 5% of
  the reads) reached the report as A297M from FDSTOOLS. Now FLASH has R1's
  base there as well and keeps the one with the higher quality; the error is
  under 1% in every sample. R1 still stops before the poly-C tract.
- In the merge, rCRS on fewer than `--minimum` reads in Mutect2's counts no
  longer counts. Where the two callers' substitution labels differ, the merge
  averages their rCRS shares to decide whether rCRS is still present, and
  takes Mutect2's share from its read counts, where a single rCRS read
  counted although Mutect2's own processing ignores it: FDSTOOLS with rCRS
  on 2.5% and Mutect2 with 1 rCRS read of 21 averaged to 3.6% and gave
  T16189Y. Now that read counts as 0, as everywhere else.
- Mutect2 no longer reports C3106c, a deletion no molecule carries. 3107 is
  rCRS's placeholder N, a base no molecule has. A read can be aligned with its
  gap at 3106 instead and its C on the N; Mutect2 then records the C at 3106
  as deleted, and once a few reads are aligned that way the report showed
  C3106c at about 50%. That deletion is now dropped from Mutect2's records,
  as FDSTOOLS' N3107DEL is; no mitoLEAF haplogroup deletes 3106.
- A sample whose Mutect2 records are all dropped (a negative control whose
  only record is a length variant under the length-heteroplasmy floor or on
  too few reads, for one) gets an empty Mutect2 table. The step stopped with
  an error and halted the run.
- Mutect2's insertions and deletions are placed using only the substitutions
  the whole sample carries. To place an insertion or deletion at the end of
  its run, the post-processing writes the sample's substitutions into a copy
  of rCRS. It wrote every substitution, even one on a single read: one read
  with T310C then joined the C-runs 303-309 and 311-315, and the C inserted
  after 302 in the other molecules was placed at 315.1C instead of 309.1C (and
  added to the real 315.1C, above 100%).
- Mutect2's length variants are added up before the length-heteroplasmy floor
  and the read minimum apply. Mutect2 writes one record per length of a run
  (one, two or three extra Cs), and a report label such as 309.1C is the sum
  of every record that carries it. Records under the floor or on fewer than
  `--minimum` reads were dropped one at a time, before adding up: a record
  with three extra Cs on 5.5% was missing from 309.1C (94.3% instead of
  99.8%) and from 309.2C. Now each label's total is checked, as for FDSTOOLS.
- A Mutect2 deletion whose records add up to the length-heteroplasmy ceiling
  (90%) is major (A523-), as an insertion already was (309.1C). One AC unit
  deleted on 85% of the reads and two on 10% is a deletion of 523-524 on 95%;
  the report wrote A523a C524c (minor). The Type column now says DEL or INS
  for a major length variant; it kept the Type of the first record (LHP).
- Mutect2's deletions are placed by the sample's own bases. Mutect2 writes
  the rCRS bases a deletion removes, and the post-processing took those
  letters to place the deletion at the end of its run. Where the sample has
  its own base in the deleted stretch, that is wrong: with T16189C, TC deleted
  after 16188 removes two Cs from a run of ten, so the deletion is C16192c
  C16193c. A rule for 16188 alone counted one base instead of two (C16193c),
  and elsewhere rCRS's letters were placed: with T310C, two Cs deleted after
  309 came out as T308t C309c instead of C314c C315c.
- Mutect2's records are no longer filtered in p12. p12 kept a record only if
  its frequency reached `--min_vf_MT2`, one record at a time, so the alleles
  of a length variant were filtered before they were added up: an allele with
  two extra Cs in the C-run 303-309 on 4.7% was missing from 309.1C.
  `--min_vf_MT2` now applies in the post-processing, after adding up, as for
  FDSTOOLS.
- The primer BED (for IGV) matches the primer files and the reads. Amplicon
  045's left primer starts at 7086, not 7087; 046's second left oligo
  (7272-7297, T at 7274) is added; 100's right primer is drawn at 34-53 on
  the appended copy (16603-16622), where its reads align, instead of at 1-20.
- The library's flank length for amplicon 001's left flank is 18, the length
  of the flank (it said 18-23; FDSTOOLS used 18 either way). A comment says
  why the flank, like the primer file, begins at position 1 although the
  oligo begins 5 bases earlier.
- FDSTOOLS reads amplicon 100 up to position 33, where its right primer
  begins, with the two primers as flanks, as for every other amplicon. The
  library ended 100 at 31 and had its flanks commented out, so FDSTOOLS took
  16 reference bases on each side as flanks and counted 32-33 as flank.
- Amplicon 045's left primer ends at 7112. The primer file gave it as 46 bases
  (7086-7131), so `p05` cut 7113-7130 from every 045 read along with the
  primer, and Mutect2 saw those positions only in amplicon 044. The primer of
  the earlier chemistry (7107-7130, `045_pA_7132-7310_l_short`) stays, so data
  made with it are trimmed as before. FDSTOOLS keeps 045 from 7131, where that
  primer ends, so it can read both chemistries; a comment in the library says
  so.
- cutadapt (`p05`) accepts a primer that is incomplete at the end of a read.
  Where the first bases of a primer's oligo carry an error, bwa clips them and
  `p03` removes them, so the read begins a few bases inside its primer; and
  quality trimming ran before the primer search and could cut into the right
  primer. cutadapt required the whole primer at each end and discarded these
  reads. The primers are now non-internal (`X` instead of `^` and `$`, set in
  `p05`; the primer files are unchanged): a primer may be incomplete at the
  read end, with at least `params.primer_min_overlap` (15) of its bases in
  the read, but is never searched inside the read. Quality trimming runs
  after the primers are removed. Mutect2 gets about 2% more reads and
  FDSTOOLS about 1%; no call changes.
- Read-through past the primers is removed in `p03` even where it matches the
  reference. A read that runs past the end of its amplicon continues into the
  kit's own sequence; `p03` removed that only where it did not align. Where the
  first bases happened to equal the mtDNA next to the primer, they stayed, and
  merged reads began up to a few bases before the left primer: cutadapt then
  did not find the primer and discarded the read (16% of amplicon 049's
  reads, 5% of 075's). `p03` now cuts a read back to the nearest primer
  boundary it passes by at most 5 bases; the boundaries come from the primer
  files cutadapt uses. A product spanning two amplicons begins and ends at
  primer boundaries and stays whole.
- `rCRS_NimaGen.fasta` trimmed from 16623 to 16622 bp, so the appended copy of
  chrM:1-53 ends where the origin-spanning amplicon's primer does.
- The humans-reference BWA index is built once, in a new
  `p01b_prepare_humans_index` process, instead of inside `p02` for every
  sample, where parallel samples could race to build it on a first run.
  `bunzip2 -k` now keeps the source archive. It also builds the humans
  file's `.fai`, which rtn otherwise creates itself, so parallel rtn runs on a
  first run could race to write it.
- `params.mapQ` now applies to the BAM Mutect2 runs on: `p09` built the
  MAPQ-filtered BAM but passed the unfiltered one on. Mutect2 already ignores
  reads below MAPQ 20 by default, so calls are largely unaffected; QC and the
  published BAM now match what Mutect2 sees. The unfiltered BAM is kept as
  `*.rtn.unfiltered.bam`.
- Left primers for amplicon 045 include its longer primer variant, so reads
  amplified with it are no longer discarded as untrimmed.
- Mutect2 deletions in repeats were shifted one base short of their rightmost
  position (`shift_deletion_right`) and could disagree with FDSTOOLS by one
  position (e.g. chrM:8281-8289 reported as 8280-8288).
- Mutect2 variants seen through the appended chrM:1-53 copy are labelled with
  true rCRS positions (`T10C`, not `T16579C`). Previously only 16570-16587 was
  mapped back, and only for sorting. Where both amplicons cover the same bases
  (chrM:21-33), the two Mutect2 records are merged by summing their reads.
- Minor deletions are labelled the same way by both callers (lowercase ref
  base, e.g. `A523a`). FDSTOOLS wrote every deletion as `A523-` regardless of
  frequency, so one deletion could appear on two rows.
- The merge step puts a variant spelled differently by the two callers on one
  row instead of two (`-309.1C`/`-309.1c`, `T16189C`/`T16189Y`,
  `A523-`/`A523a`).
- `LOW` (FDSTOOLS) now flags an amplicon when its reads, summed over all its
  rows except "Other sequences" (the reads frequencies are computed on), are
  below `--depth`, checked for every amplicon in the library. Previously only
  amplicons reported as a single row were checked, so one with a sequence
  plus "Other sequences" was never flagged.
- Different Mutect2 alleles at one position are reported as separate rows.
  They were merged into one row with frequencies and read counts added and
  only the first label kept (C>A 29.3% and C>T 5.1% became `C756M` at 34.4%).
  Identical claims, such as two alleles of one insertion, are still combined.
- FDSTOOLS frequencies are divided by the amplicon's counted reads (without
  "Other sequences"). The depth was reconstructed from each row's rounded
  percentage, which slightly overestimated it and lowered frequencies; at low
  depth a call on every read could show as 95% instead of 100%.
- rtn runs with `-i`: indels no longer count toward a read's distance from
  the human genomes. rtn chose among equally close genomes by substitutions
  only but then judged the read with indels counted, so a read was kept or
  dropped depending on the order bwa listed those genomes. Reads carrying a
  length variant that no genome has are no longer dropped for it, which
  raises minor length variants in the C-stretches slightly.
- FDSTOOLS is pinned to 2.1.1 in `FMP-NimaGen.yml` and the Dockerfile. In
  2.2, `samplestats` reads every allele name back as a sequence and rejects
  `N3107DEL`, which FDSTOOLS itself writes for the rCRS placeholder at 3107,
  so `p11` failed on a fresh install.
- Building the humans index no longer makes rtn rerun. `p01b` writes the index
  into `--humans_index_dir`, which changes the folder's time stamp; on the
  next run `p01b` ran again, in a new work folder, and `p08`/`p09`, which got
  the folder through it, ran again for every sample. They now take the folder
  directly and only the file name from `p01b`, once the index is ready.
- Two runs on the same reads give the same reports. FLASH ran with one thread
  per processor and wrote the merged reads in a different order on every run;
  a few reads at the limit of the NUMT filter were then kept in one run and
  dropped in another, so read depths could differ by a few reads. `p04` now
  runs FLASH with `-t 1`.
- R1 is cut before the poly-C tract at 303-315 (`p02`) only where the whole
  anchor sequence (chrM:278-302) is found. With cutadapt's default minimum
  overlap of 3, any R1 read that happened to end like the start of the anchor
  also lost its last 3 or more bases, and the reads of the amplicon that ends
  at 296 lost their last 19.
- `tests/simulation/score_calls.py` checks all mismatches within one C-stretch
  together: A16182C against A16182- -16193.1C (one molecule, two spellings)
  lie 11 bases apart and were checked separately, as a wrong sequence. Its
  summary gives the samples whose calls describe the true sequence.
- Every `publishDir` overwrites: with `-resume`, Nextflow left a file already
  in `--outdir` in place, so after a change of option (e.g. `--frame shared`,
  then back to `--frame separate`) the folder kept the report of an earlier
  run although the step had run with the new setting or come from the cache.
- README: the repository is `FMP-NimaGen` (was `mtDNA-NimaGen`), the conda
  environment is activated as `FMP-NimaGen` (not `FMP-NimaGen.yml`), the
  overview lists the steps in the order the pipeline runs them, with the NUMT
  filter, Conda is the recommended way to run it (the Docker image predates
  this version and is untested with it), and Nextflow's cache folder is
  `.nextflow` (not `.nextflow.cache`).
- The other tools in `FMP-NimaGen.yml` and the Dockerfile are pinned too
  (gatk4 4.6.2.0, bcftools 1.21, pandas 2.2.3, biopython 1.85, matplotlib
  3.9.4, openpyxl 3.1.5, pytest 8.4.2), the versions the pipeline is
  validated with; a fresh install took whatever was newest. pytest comes
  from pip: conda's 8.4.2 needs Python 3.10.
- Outside the frame regions, FDSTOOLS no longer writes a substitution and a
  deletion at one position as one lowercase label with their shares added:
  A on 40% and 30% deleted at 756 gave C756a 70%, A on 70% and 30% deleted
  C756a 100%. They are two rows, as in the frames: the bases present (C756M
  40%, or C756A 70%) and the deletion (C756c 30%).

### Changed
- The parameter log (`p00_parameters.txt`) has shorter lines: paths inside
  the pipeline folder are written relative to it, and the longest
  explanations are shorter. It also lists the anchor before the poly-C tract
  (`--adapter`), the primer files and the amplicon positions file.
- In 16180-16193 and 300-315, one base inside an A-run that still has all its
  A's is written as an insertion, in both frames. AAACA at 16180 (four A's
  with a C between the third and the fourth) was A16183C C16184A and is now
  -16182.1C; the rest of the molecule is read as a full A-run and its C's.
  FDSTOOLS and Mutect2 align the C as an insertion, but Mutect2's main calls
  were rewritten in the old spelling: where the C-run varies in length, its
  length variants matched neither its main calls nor FDSTOOLS' rows
  (-16193.1C at 99.9% from Mutect2, -16193.1c at 34% from FDSTOOLS, for the
  same molecules). Both callers now give the same rows. M73a, M73a1 and A2g2
  are written the same way (AACAA: -16181.1C C16188-), so for them the shared
  frame no longer matches mitoLEAF, which writes 16184A. An A-run that has
  lost one of its own A's (AACAC) stays A16182C.
- Amplicons with fewer than 20 reads are marked LOW by default (`--depth`;
  before 10), the recommended minimum read depth. The calls in such an
  amplicon stay in the report; LOW only warns that they rest on few reads.
- A base counts as present from 2% of the reads by default (`--min_vf_FDS`
  and `--min_vf_MT2`; before 5%), so that low-level point heteroplasmies are
  detected. A base still needs at least `--minimum` reads (2), so where an
  amplicon has few reads, a single read cannot make a call.
- Every base in a Mutect2 call, rCRS included, and every Mutect2 length
  variant needs at least `--minimum` reads (default 2), as FDSTOOLS needs for
  every sequence it reports. A base on fewer reads does not count: a minor
  base on one read is no call, and a single rCRS read no longer turns a
  homoplasmic call into a heteroplasmy.
- rtn runs once, on the primer-trimmed reads (`p09`); `p08` keeps the reads of
  the FDSTOOLS branch whose names `p09` kept. rtn judged the FDSTOOLS reads
  with their primers, which are the kit's oligos: where a primer carries a
  rare rCRS base (750A, 8860A, 14766C) every read differed from almost every
  genome in the humans file, and a read that bwa aligned only in part was
  judged on that part. The two rtn runs could decide differently on the same
  molecule, so FDSTOOLS could miss a molecule Mutect2 saw (e.g. a length
  variant in 16180-16193). FDSTOOLS now gets the reads Mutect2 gets, with
  their primers. It no longer counts reads that cutadapt rejects (no complete
  primer, or shorter than 60 bases), and assigns about 0.7% fewer reads to
  amplicons.
- cutadapt's maximum length (`params.maximum_length`) is 600 (was 300),
  longer than any merged read from 2x226 sequencing. A product spanning two
  neighbouring amplicons of one primer pool is up to about 390 bases after
  primer trimming and was discarded before Mutect2. It is a real product of the
  mtDNA, aligns end to end, and now adds up to a quarter more reads to
  amplicons 001 to 003 and 5-8% to 098 to 100. FDSTOOLS does not use this
  limit.
- The report column `interp_total` is now `rd_total_FDS`: all reads FDSTOOLS
  counted at the position, the number `rd_FDS` is a part of and `vf_FDS` is
  calculated from. Nothing is interpolated any more; only the name changes.
- `main.nf` runs on Nextflow 26.04, which reads scripts with its strict syntax
  and refused the file: the parameter log and the three input checks moved
  from the top of the file into the workflow block, `assert` takes its message
  after a colon, a `publishDir` path that names the sample is a closure, and
  channel factories are called through `channel`. `nextflow.config` names the
  versions the pipeline is tested with (25.04.6 and 26.04.6). Output is
  unchanged.
- Python 3.12 (was 3.9, which no longer gets security fixes), samtools and
  bcftools 1.24 (were 1.21), bwa 0.7.19 (was 0.7.18) and cutadapt 5.2 (was
  4.6), in `FMP-NimaGen.yml` and the Dockerfile. The reports are the same on
  the real and simulated samples tested; FDSTOOLS counts a read or two more
  in a few amplicons.
- The report is named after its frame, `<sample>_separate_frame.xlsx` or
  `<sample>_shared_frame.xlsx` (was `<sample>_merged_variants.xlsx`), so
  reports in both frames can sit side by side; `variant_note` no longer says
  "separate frame" or "shared frame" on every row of the C-stretches.
  `tests/simulation/score_calls.py` takes `--frame`.
- rtn compares reads with `humans_NimaGen_2026-10.fa` (built by
  `resources/rtn_files/humans/build/` from public sources) instead of rtn's
  `humans.fa`: the mitoLEAF haplogroups, the GenBank genomes on mitoLEAF's
  list and the re-called 1000 Genomes and HGDP people, without genomes that
  carry NUMT-like stretches or errors. Fewer reads of rare lineages are taken
  for NUMTs, and no more NUMT reads pass. `p01b` unpacks the archive named by
  `--humans_index_base` (it always unpacked `humans.fa.bz2`), and only if it
  is not unpacked yet; the new name makes the first run after the update build
  the new index once.
- Mutect2's major calls in 57-60, 300-315 and 16180-16193 and within 12 bases
  of them are written as FDSTOOLS' molecules are: rebuilt into one molecule,
  placed by the same alignment, and written in the frame chosen with `--frame`
  (edges by the general rule), so both callers meet on the same labels.
  Mutect2 wrote D5c's molecule as `A16181- A16182- A16183- -16192.1T ...` and
  T57C with one T more as `-56.1C`. Minor calls there stay as Mutect2 writes
  them.
- Major calls outside 57-60, 300-315 and 16180-16193 are written by the
  general rule (`resources/scripts/notation.py`) in both callers' tables:
  calls within 10 bases of each other are applied to rCRS and described anew.
  One molecule could otherwise be spelled two ways and land on two rows
  (`T15940- T15941C` from one caller, `-15939.1C T15943- T15944-` from the
  other, both `T15940C T15944-`). A base at the rCRS N is now `-3109.1T`
  instead of `N3107T`. Minor calls stay as each caller writes them; a label the
  rewrite keeps keeps its frequency, new labels take the lowest frequency of
  the rows they replace, and rewritten FDSTOOLS rows name the original labels
  in `variant_note`.
- In 57-60, 300-315 and 16180-16193, FDSTOOLS' rows come from the frames:
  every sequence of the amplicon is placed into the region and written in the
  frame chosen with `--frame` (`separate` by default, or `shared`; 57-60 is
  laid from the left in both). This replaces FDSTOOLS' own labels there, which
  could merge a substitution and a deletion into one label (`T16189c`) and
  wrote T57C with one T more as `-56.1C`. The rows carry the frame in
  `variant_note`. Within 12 bases of each region the labels come from the same
  alignment, so a change next to the region is not counted on both sides of
  its edge.
- Mutect2 always evaluates the alleles in `--force_alleles`
  (`resources/mutect2/force_alleles.vcf`, only A3105G for now). A3105G sits
  two bases before the rCRS N at 3107, where bwa places the missing base as a
  deletion plus a substitution; with pileup detection alone Mutect2 reported
  nothing there. In samples without it the forced allele comes out near 0% and
  is filtered. Forcing alleles inside dense clusters of variants made Mutect2
  lose calls there, so the list is kept to isolated sites.
- Mutect2 also takes candidate changes straight from the aligned reads
  (pileup detection, substitutions and indels seen in at least 10% of the
  reads), not only from its local assembly, which missed variants in dense
  clusters, next to C-stretches and close to amplicon ends. Only bases of at
  least `--baseQ` count (Mutect2's default is 12), and the proper-pair check is
  off because merged reads are unpaired and would all be rejected.
- Mutect2 reports each substituted position as its own record
  (`--max-mnp-distance 0`). By default it merged substitutions on neighbouring
  positions carried by the same reads into one multi-base record, which split
  each position's frequency between records and reached the report as one row
  with one frequency for several positions.
- `p07` maps the primer-trimmed merged reads with `bwa mem -L 100,100`
  (`--clipping_penalty`; bwa default 5,5). After primer trimming, a variant a
  few bases from an amplicon end sits at the end of every read covering it,
  and the default clipping cut it off, so Mutect2 missed it. `p02` keeps the
  default: there clipping removes adapter and off-target read ends for `p03`.
- Length heteroplasmy uses one symmetric threshold (`--lh_thresh`, default
  10): below 10% not reported, 10-90% reported as minor (lowercase), 90% and
  above as major. Previously there was no lower floor and the default of 90
  only set the major cut-off; 10 and 90 now give the same result. Applies to
  insertions and deletions in both callers.
- When both callers report a variant but disagree on major vs minor, the
  reported call follows their averaged frequency, and the caller on the other
  side is flagged `DISAGREEMENT` (own colour in the Excel output). `False` now
  only means the caller did not report the variant.
- `LOW` rows cover both callers: one row per amplicon, at the bottom of the
  merged table. The FDSTOOLS and MUTECT2 columns name the amplicon when that
  caller's depth is below `--depth`, and `called_by_*` read `LOW` (red, like a
  missed call) or `OK` for each caller, instead of True/False, which in variant
  rows mean called or missed. Mutect2's depth is the read count at the
  amplicon's middle position in the BAM it runs on (from `p09`).
- `--min_reads_per_strand` default 3 → 0. Merged reads are single-orientation,
  so any non-zero value tagged every Mutect2 call `strict_strand`. The FILTER
  column is informational; no calls were dropped by it.
- The IUPAC codes live in one table, `resources/scripts/iupac.py`, shared by
  both callers' scripts, the frames and the merge (each kept its own copy of
  the two-base codes). It also holds the three- and four-base codes; output
  is unchanged.
- The merge compares substitutions by position: each caller's call at a
  position is one row, whatever bases it found (T16189C, T16189Y, T16189H).
  Where the callers' labels differ, every base other than rCRS that either
  found is kept, and rCRS counts as present when the callers' averaged share
  of it (molecules with a deletion there left out) reaches `--min_vf`; a
  caller whose label differs is DISAGREEMENT. Different bases from the two
  callers were two rows, each missed by the other caller. A caller's rows
  per base at one position are read as one row (T16189Y and T16189W become
  T16189H). With several bases other than rCRS, vf and reads are given per
  base (`C 40, A 10`); such values no longer turn the column into text, which
  broke the averaging on the other rows. Three- and four-base codes are
  coloured like the two-base ones.
- FDSTOOLS writes one row per position for substitutions, from its labels and
  from the frames' rows: the bases with at least `--min_vf_FDS` of the reads
  are present, rCRS when the reads without another base or a deletion there
  reach it. One base other than rCRS present is a major call (C756A), two or
  more give their IUPAC code, three and four bases included (T16189H for T, C
  and A; T16189M for C and A without T), with vf and reads per base other
  than rCRS (`C 40, A 10`). Several bases were one row each, coded with rCRS
  even where rCRS was absent, and a base was major only from 95%: C on 94%
  of the reads with T and A on 3% each was T16189Y, now T16189C.
- Mutect2's substitutions are one row per position by the same rule as
  FDSTOOLS', from `--min_vf_MT2`: its records for several bases at a position
  (A and T at 756: C756M and C756Y) are one row (C756H, vf `A 0.293, T 0.051`,
  reads `70,29,5`). Whether rCRS is present is decided from Mutect2's reads
  (its reference reads among all reads it counted there), here and in the
  merge, not from 1 minus its frequency: for a base on all n reads Mutect2
  reports about (n + 1) / (n + 2), so with `--min_vf_MT2 2` every call from
  fewer than about 48 reads kept rCRS and became a code (with 1%, 100 reads
  all G gave A73R instead of A73G). The frequencies shown stay Mutect2's
  (README, Known Limitations).

### Added
- README, Reading the Report: Mutect2 also uses reads that span two
  neighbouring amplicons of one primer pool, which FDSTOOLS cannot assign, so
  its depth is higher there.
- `--variant_note false` leaves the `variant_note` column out of the report;
  by default it stays.
- README, Known Limitations: recurring low-quality sequencing errors can
  appear as rows from FDSTOOLS only, which does not use base qualities; a
  higher detection threshold leaves them out.
- README, Known Limitations: where Mutect2's rows in 16180-16193 differ from
  FDSTOOLS' (an A-run one A short with C-runs of different lengths; leave
  Mutect2 out there with `--mutect2_disabled_regions`), and that `vf_MT2` is
  Mutect2's own estimate, not the share of its reads.
- `resources/primers/README.md`: which file describes the primers for which
  step, how their ends were checked against the reads, what happens to reads
  from a primer bound at a wrong site, what the wobble codes mean, and why
  amplicons 001, 045, 046 and 100 need a note.
- README: Install rtn, from the fork at <https://github.com/PeterResutik/RtN>
  (branch `fixes`: base qualities used, a tie between hits broken by the
  likelihood rtn decides on, inserted bases counted as aligned with `-i`) with
  its `build.sh`.
- `--rtn_positions` (by default the 50 positions in
  `resources/rtn_files/positions/`, `none` to switch it off): `p09` gives rtn a
  copy of the reads with base quality 2 at these positions, keeps rtn's decision
  and gives the reads their own qualities back, so Mutect2 sees the original
  ones (`resources/scripts/rtn_qualities.py`). With an rtn that uses base
  qualities, a difference at these positions counts about a tenth of a full
  mismatch; rtn as released ignores qualities, and nothing changes.
- `resources/rtn_files/positions/`: the 50 positions that change most often in
  human mtDNA (16519, 16182 and 16183, then mitoLEAF 1.6.1's tree by number of
  branches with a substitution), for rtn to count at base quality 2. The number
  was chosen by a rule fixed in advance on a leave-out of the humans file
  (`curve.tsv`); the README gives the rule, the results and how to rebuild the
  list with `rank_positions.py`.
- `p00_parameters.txt` names the human mitogenomes and NUMT sequences rtn
  compares reads with (`--humans_index_dir`, `--humans_index_base`,
  `--numts_index_dir`, `--numts_index_base`), so every run records which
  humans file it used.
- README: Input and Output (file names, the results folders) and Main Options.
- README: Reading the Report (notation, columns, flags and colours).
- `resources/amplicon_bed/NimaGen_primers.bed`: the panel's 200 primer binding
  sites and 64 wobble bases, for checking primer sites against population
  variants. The pipeline does not read it.
- `dominant_molecule` column on the rows of 16180-16193, 300-315 and 57-60:
  the change of the region's dominant molecule at the row's position and of
  its kind (substitution, deletion, or insertion at that index), in the
  selected frame, with the molecule's share of the reads, e.g. `-16193.2C
  (40.9%)`. The dominant molecule is FDSTOOLS' molecule with the most reads;
  on equal reads the one with fewer changes from rCRS; molecules equal on
  both are shown together (`T16189C (25.0% + 25.0%)`). Empty where the
  dominant molecule does not carry that change. The rows still give every
  molecule's share per position; the column shows which changes belong to
  the most common molecule.
- `LICENSE`: MIT, copyright University of Zurich. Third-party material in the
  repository keeps its own terms, listed in the README.
- `--skip_numt_filter` (default `false`): skips the NUMT filter (`rtn`) in
  `p08`/`p09`, for simulated reads, which contain no NUMTs. The MAPQ filter
  and everything downstream are unchanged.
- `--publish_bams` (default `true`): with `false`, BAM/FASTQ intermediates of
  `p02`, `p03`, `p08`, `p09` and Mutect2's bamout are not copied to the output
  folder; VCFs, FDSTOOLS files, QC and the merged Excel still are.
- `tests/simulation/`: tools that simulate NimaGen reads for mitoLEAF
  haplogroups, run the pipeline on them and score the calls against the
  known sequence (see its README).
- `simulate_reads.py --ambiguous alt` and `rcrs`: positions mitoLEAF marks as
  ambiguous always take the variant or the rCRS state, besides one state drawn
  per haplogroup (`one`, the default) and a 50/50 heteroplasmy (`het`).
- `resources/scripts/frames.py`: the shared and the separate frame for
  16180-16193 and 300-315 (57-60 is laid from the left, with no separate
  variant), and the report rows they give per sample. FDSTOOLS' molecules are
  placed into a region by aligning the amplicon's whole sequence to rCRS, so a
  variant in the flanks cannot break the cut. Not used by the report yet.
- `resources/scripts/notation.py`: how a molecule differs from rCRS, by the
  general rule (fewest changes, substitutions before gaps, gaps placed 3',
  whole repeat copies at the 3'-most copy, the rCRS N at 3107 left out). Not
  used by the report yet; the frames for the complex regions and the rebuild
  of the callers' major calls will build on it.
- Unit tests for the Python report scripts in `tests/`, run with `pytest` from
  the repository root (`pytest` added to `FMP-NimaGen.yml`). The first file
  covers `process_mutect2_output_improved.py`: origin coordinates and pooling,
  the length-heteroplasmy floor and ceiling, minor labels, deletions in
  repeats, and different alleles at one position kept apart.
- Tests for `process_fdstools_output_improved_better.py`: frequencies counted
  on an amplicon's reads without "Other sequences", LOW for low or missing
  amplicons with their calls kept, the length-heteroplasmy floor and ceiling,
  and minor deletions and substitutions.
- A test that the frame labels of every run structure of 16180-16193 and 300-315
  (run lengths, interrupt present, missing or substituted) and of 57-60 give
  their molecule back exactly, in both frames.
- Tests for `merge_fdstools_mutect2_improved.py`: one row per locus across
  both callers' spellings, DISAGREEMENT kept apart from a missed call, the plain
  and depth-weighted averages, LOW rows per amplicon for both callers, and the
  colours of these flags in the Excel file.
- `resources/rtn_files/humans/build/`: the scripts, lists and calls table
  that build `humans_NimaGen_2026-10.fa`, rtn's human mitogenomes for NimaGen
  (`build.sh`; sources, versions, screening rules and licences in its
  README). Each genome is written once, plus its first 100 bases.
- `--disagreement_average` (`plain` by default, or `depth_weighted`): how the
  two callers' frequencies are averaged when they disagree on major vs minor.
  `depth_weighted` weights each caller by its read depth for the call, so a
  call on a few Mutect2 reads no longer outweighs one on many FDSTOOLS reads.
- `--mutect2_disabled_regions` (`none` by default, `all`, or region names
  separated by commas, e.g. `300-315,16180-16193`): complex regions where only
  Mutect2's major calls count, as rebuilt in the frame. Its minor calls there
  (PHP, LHP) leave the report and are listed in the `variant_note` of the
  region's rows (not shown when FDSTOOLS reports nothing in the region); a row
  there that Mutect2 has no call for says `DISABLED` (grey) instead of False. Applied
  in the merge, so with `-resume` only `p13` reruns.

### Removed
- `--initial_tumor_lod`, `--tumor_lod_to_emit`, `--max_reads_per_alignment_start`
  and `--min_reads_per_strand`: their value (0) is fixed in `p12`, since other
  values skew or block the calls for these amplicon reads. The parameter log
  lists them with Mutect2's other fixed options.
- rtn's `humans.fa.bz2`: unused since rtn compares reads with
  `humans_NimaGen_2026-10.fa`; it stays in the repository's history.
- `--detection_limit`: unused since the mutserve-based setup was replaced; the
  Mutect2 frequency floor is `--min_vf_MT2`.
- `--min_reads_filt`: only passed to `fdstools samplestats`, where it has no
  effect because samplestats filtering is off.
- `csvtk` from `FMP-NimaGen.yml` and the Dockerfile: no step of the pipeline
  uses it.
- FastQC from `p10` and the environment: nothing used its report, and on
  amplicon reads it flags duplication and sequence content on every sample.
  `p10` keeps the read-depth plot.

## [0.1.4] – 2025-05-19

### Fixed
- Added `openpyxl` to the Conda environment to enable Excel output.
- Resolved Apple Silicon compatibility issue with `cutadapt` by falling back to `pip install cutadapt` in environments where the Conda build fails.

## [0.1.3] – 2025-05-16

### Fixed
- Added missing tools `cutadapt=4.6` and `fdstools` to `FMP-NimaGen.yml`
- Ensured pipeline now runs end-to-end without missing dependencies

## [0.1.2] – 2025-05-16

### Fixed
- Corrected the filename of the Conda environment file (`FMP-NimaGen.yml`)
- Updated usage instructions in `README.md` to match filename
- Changed output directory path in `main.nf` from `results_new_7/` to `results/`

## [0.1.1] – 2025-05-16

### Fixed
- Removed `.ipynb_checkpoints` folders from repository
- Added `.ipynb_checkpoints/` to `.gitignore`
- Clarified that Jupyter notebooks were only used for development


## [0.1.0] – 2025-05-16

### Added
- First public release of the **FMP-NimaGen** pipeline.
- Support for forensic mtDNA variant calling using **FDSTools** and **Mutect2**.
- Integrated merging of variants from both callers, including caller-specific flags.
- Variant formatting compliant with **EMPOP** standards.
- Resolution of point and length heteroplasmies using IUPAC and lowercase notation.
- Position normalization for circular mtDNA: **16570–16587 → 1–18**.
- Conda environment specification for reproducible local runs.
- Docker image for Intel-based systems (Linux, macOS, and Windows with x86_64 architecture).


### Known limitations
- Pipeline currently supports **only NimaGen mtDNA panel data**.
- **Docker not yet compatible with Apple Silicon (ARM64)**.
- No multi-panel or profile support (e.g., Thermo Fisher or custom kits).

