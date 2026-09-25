# Changelog

All notable changes to this project will be documented in this file.

This project follows [Semantic Versioning](https://semver.org) and uses the [Keep a Changelog](https://keepachangelog.com) format.

---

## [Unreleased]

### Fixed
- **Race condition in humans-reference BWA index build**: the check-and-build
  logic for the humans reference index lived inline in `p02`, which runs once
  per sample — on a first run, parallel samples could all see the index as
  missing at once and race to `bunzip2`/`bwa index` the same file
  concurrently. Moved the logic into a new one-shot process
  (`p01b_prepare_humans_index`) that `p02`, `p08`, and `p09` now depend on,
  so it runs exactly once regardless of sample count. Also changed
  `bunzip2` to keep the source `.bz2` (`-k`) instead of deleting it.
- **MAPQ filter never reached Mutect2**: `p09` computed a MAPQ≥`params.mapQ`
  (default 30) filtered BAM (`samtools view -q ...`) but only used it to
  generate a QC read-depth statistic, then discarded it — the BAM actually
  propagated to Mutect2 (`p12`) was the *unfiltered* NUMT-filtered output of
  `rtn`. This let low-confidence/ambiguous reads (e.g. MAPQ 0) reach the
  variant caller uncontrolled, which can masquerade as low-frequency
  (heteroplasmic) variant signal. It also meant `p10`'s FastQC ran on a
  different BAM than the one its read-depth plot was computed from. `p09`
  now indexes and propagates the mapQ-filtered BAM as its `bam_file` output;
  the unfiltered RTN output is kept alongside (renamed `*.rtn.unfiltered.bam`)
  for comparison/debugging.

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

