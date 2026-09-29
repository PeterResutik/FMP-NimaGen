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
  `bunzip2 -k` now keeps the source archive.
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

### Changed
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
  caller's depth is below `--depth`, and `called_by_*` say which caller is
  low. Mutect2's depth is the read count at the amplicon's middle position in
  the BAM it runs on (from `p09`).
- `--min_reads_per_strand` default 3 → 0. Merged reads are single-orientation,
  so any non-zero value tagged every Mutect2 call `strict_strand`. The FILTER
  column is informational; no calls were dropped by it.

### Added
- `--disagreement_average` (`plain` by default, or `depth_weighted`): how the
  two callers' frequencies are averaged when they disagree on major vs minor.
  `depth_weighted` weights each caller by its read depth for the call, so a
  call on a few Mutect2 reads no longer outweighs one on many FDSTOOLS reads.

### Removed
- `--detection_limit`: unused since the mutserve-based setup was replaced; the
  Mutect2 frequency floor is `--min_vf_MT2`.
- `--min_reads_filt`: only passed to `fdstools samplestats`, where it has no
  effect because samplestats filtering is off.

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

