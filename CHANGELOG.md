# Changelog

All notable changes to this project will be documented in this file.

This project follows [Semantic Versioning](https://semver.org) and uses the [Keep a Changelog](https://keepachangelog.com) format.

---

## [Unreleased]

### Fixed
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
- The report is named after its frame, `<sample>_separate_frame.xlsx` or
  `<sample>_shared_frame.xlsx` (was `<sample>_merged_variants.xlsx`), so
  reports in both frames can sit side by side; `variant_note` no longer says
  "separate frame" or "shared frame" on every row of the C-stretches.
  `tests/simulation/score_calls.py` takes `--frame`.
- rtn compares reads with `humans_NimaGen.fa` (built by
  `resources/rtn_files/humans/build/`) instead of rtn's `humans.fa`. NUMT
  reads that one of the removed genomes let through in the control region are
  now filtered, reads of rare lineages are no longer dropped in amplicons no
  genome matched, and the file and its bwa index are half the size. `p01b`
  unpacks the archive named by `--humans_index_base` (it always unpacked
  `humans.fa.bz2`), and only if it is not unpacked yet; the first run after
  the update builds the new index once.
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
- Mutect2 records resting on fewer than `--depth` reads (FORMAT/DP) are
  dropped. In the merge, a Mutect2 frequency from one or two reads counted as
  much as an FDSTOOLS call from many (a 1-read 67% turned a 22-read 100%
  `T16189C` into `T16189Y`). The amplicon's LOW row still flags the low
  coverage.
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
- README: Input and Output (file names, the results folders) and Main Options.
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
- `resources/rtn_files/humans/build/`: scripts and lists that build
  `humans_NimaGen.fa`, rtn's human mitogenomes for NimaGen, from rtn's
  `humans.fa`: four genomes left out (two containing NUMT sequence, two not
  modern human), 580 mitoLEAF haplogroups added where no genome matched a
  NimaGen amplicon exactly, and each genome written once plus its first 100
  bases instead of twice.
- `--disagreement_average` (`plain` by default, or `depth_weighted`): how the
  two callers' frequencies are averaged when they disagree on major vs minor.
  `depth_weighted` weights each caller by its read depth for the call, so a
  call on a few Mutect2 reads no longer outweighs one on many FDSTOOLS reads.
- `--mutect2_disabled_regions` (`none` by default, `all`, or region names
  separated by commas, e.g. `300-315,16180-16193`): complex regions where only
  Mutect2's major calls count, as rebuilt in the frame. Its minor calls there
  (PHP, LHP) leave the report and are listed in the `variant_note` of the
  region's rows (not shown when FDSTOOLS reports nothing in the region); a row
  there that Mutect2 has no call for says `DISABLED` instead of False. Applied
  in the merge, so with `-resume` only `p13` reruns.

### Removed
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

