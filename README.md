# Forensic mtDNA Pipeline (FMP)

This pipeline processes mitochondrial DNA (mtDNA) sequencing data generated using the Nimagen RC-PCR Kit. It is implemented using Nextflow, enabling streamlined execution across different computational environments.

## Table of Contents

- [Overview](#overview)  
- [Installation and Setup](#installation-and-setup)  
  - [Nextflow Requirements and Installation](#nextflow-requirements-and-installation)  
  - [Clone the Repository](#clone-the-repository)  
  - [Organize Input Data](#organize-input-data)  
- [Running the Pipeline](#running-the-pipeline)  
  - [Running Locally with Conda (Recommended)](#running-locally-with-conda-recommended)  
  - [Running with Docker](#running-with-docker)  
- [Input and Output](#input-and-output)  
- [Reading the Report](#reading-the-report)  
- [Main Options](#main-options)  
- [Configuration](#configuration)  
- [Tests](#tests)  
- [Known Limitations](#known-limitations)  
- [Cleaning Up](#cleaning-up)  
- [Citation](#citation)  
- [Contributing](#contributing)  
- [Contact](#contact)  
- [License](#license)  

## Overview

The pipeline automates the following steps:

1. **Mapping and soft-clip removal** (`p02`–`p03`): read pairs are mapped to rCRS with bwa, and soft-clipped bases are removed.  
2. **Read merging** (`p04`): FLASH merges each read pair into one read.  
3. **Primer trimming** (`p05`): Cutadapt removes the primers, trims low-quality ends and drops reads without a primer or outside the length limits.  
4. **Mapping** (`p06`–`p07`): the merged reads are mapped to rCRS with their primers (for FDSTOOLS) and without (for Mutect2).  
5. **NUMT filtering** (`p08`–`p09`): rtn removes reads that no human mitogenome in `humans_NimaGen_2026-10.fa` explains closely enough, which removes reads from nuclear copies of mtDNA (NUMTs); reads below `--mapQ` are removed too. rtn judges the reads without their primers (`p09`); `p08` keeps the same reads, by name, for FDSTOOLS. At the 50 positions that change most often in human mtDNA (`resources/rtn_files/positions/`), rtn sees the reads at base quality 2, so that differences there count little; the callers keep the original qualities.  
6. **Quality control** (`p10`): read depth per amplicon.  
7. **Variant calling** (`p11`–`p12`): **FDSTOOLS** on the reads with primers and **GATK Mutect2** on the trimmed reads.  
8. **Report** (`p13`): both callers' calls are written in one notation and merged into one Excel table per sample, which shows where the callers agree and where they do not.  

The workflow is designed for reproducibility and scalability.

## Installation and Setup

### Nextflow Requirements and Installation

Nextflow requires:

* **Bash 3.2 or later**  
* **Java 17 (or later, up to Java 24)**

Check your Java version:

```bash
java -version
````

If not installed, the easiest way is via [SDKMAN!](https://sdkman.io/):

```bash
# Install SDKMAN
curl -s https://get.sdkman.io | bash

# Restart your terminal, then install Java
sdk install java 17.0.10-tem

# Confirm installation
java -version
```

Install Nextflow:

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /usr/local/bin  # or any directory in your PATH
nextflow -version
```

### Clone the Repository

```bash
git clone https://github.com/PeterResutik/FMP-NimaGen.git
cd FMP-NimaGen
```

### Organize Input Data

Place your raw sequencing files into the expected directory:

```bash
cp /path/to/your/FASTQ/*fastq.gz raw_data/
```

## Running the Pipeline

### Running Locally with Conda (Recommended)

This is the way the pipeline is tested:

1. Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Conda](https://docs.conda.io/en/latest/) installed.
2. Create the environment:

```bash
conda env create -f FMP-NimaGen.yml
```

3. Activate the environment and run the pipeline:

```bash
conda activate FMP-NimaGen
nextflow run main.nf -profile local
```

To resume a previous run and skip already completed steps, add `-resume`.

### Running with Docker

```bash
nextflow run main.nf -profile docker
```

> **Note**: the Docker image (`peterresutik/nimagen-pipeline:latest`) dates from May 2025
> and has not been tested with this version of the pipeline. Use Conda until a new
> image is published.

## Input and Output

Read files must be named `<sample>_R1_001.fastq.gz` and `<sample>_R2_001.fastq.gz`
and are taken from `raw_data/` by default. Results go to `results/`:

| Folder | Contents |
|---|---|
| `p00_parameters.txt` | the options of the run |
| `p10_quality_control/<sample>/` | read depth per amplicon (plot and tables) |
| `p11_fdstools/<sample>/` | FDSTOOLS files (sequences per amplicon, statistics, HTML) |
| `p12_mutect2/<sample>/` | Mutect2 VCF |
| `p13_merged_variants_xlsx/<sample>_separate_frame.xlsx` | **the report** (`_shared_frame` with `--frame shared`) |

BAM and FASTQ intermediates are also copied unless `--publish_bams false` is set.

## Reading the Report

The report has one row per variant, sorted by position. The first column,
`FMP`, is the pipeline's call. The columns after it show what each caller found.

### Notation

| Example | Meaning |
|---|---|
| `A263G` | substitution: rCRS base, position, base found |
| `A523-` | deletion of the rCRS base at 523 |
| `-309.1C` | insertion of a C after 309; a second inserted base is `-309.2C` |
| `T195Y` | point heteroplasmy: two or more bases at 195, written with the IUPAC code (`Y` is C and T) |
| `-309.2c` | lowercase: a length variant on part of the molecules only (between `--lh_thresh` and 100 minus `--lh_thresh` percent) |
| `LOW` | not a variant: an amplicon with fewer reads than `--depth`, named in `marker`; no call can be made in `marker_range` |

In 16180–16193 and 300–315 the labels follow `--frame`: `separate` writes
each run of equal bases by its length, `shared` lays the sequence from the
left, as mitoLEAF and EMPOP do. Outside these regions a molecule is written
with the fewest changes from rCRS. That can differ from the spelling of the
phylogeny: L1c1c and L5b1 are written `G247-` where the phylogeny writes
`247A 249del`. The molecule is the same.

### Columns

| Column | Meaning |
|---|---|
| `FMP` | the call |
| `FDSTOOLS`, `MUTECT2` | each caller's own label; empty if the caller did not report the variant |
| `vf_FDS`, `vf_MT2` | share of the reads with the variant, in percent; with several bases at one position, one value per base (`C 40, A 10`) |
| `rd_FDS`, `rd_total_FDS` | FDSTOOLS: reads with the variant, and all reads at the position |
| `rd_MT2`, `MBQ` | Mutect2: read counts and median base quality, each for rCRS first, then for the variant (`687,1450`: 687 reads with the rCRS base, 1450 with the variant) |
| `called_by_FDSTOOLS`, `called_by_MUTECT2` | whether each caller confirms the call (see below) |
| `variant_note` | remarks, for example which Mutect2 calls were left out |
| `dominant_molecule` | in 16180–16193, 300–315 and 57–60: what the molecule with the most reads carries at this position, with that molecule's share of the reads, e.g. `-16193.1C (40.9%)`; empty if it carries no change there |
| `marker`, `marker_range`, `num_markers` | the amplicon of the FDSTOOLS call, its range, and the number of amplicons the call was seen in |
| `Filter`, `Pos`, `Ref`, `Variant`, `GT` | fields of the Mutect2 VCF record |
| `Type` | `SNP`, `PHP` (point heteroplasmy), `INS`, `DEL` or `LHP` (length heteroplasmy) |

The two callers show their evidence differently. FDSTOOLS gives the reads and
the share of the variant only (`rd_FDS`, `vf_FDS`); the reads without it
follow from `rd_total_FDS`. Mutect2 gives the reads of both the rCRS base and
the variant (`rd_MT2`), and the share of the variant (`vf_MT2`).

### Flags and colours

`called_by_FDSTOOLS` and `called_by_MUTECT2`:

| Value | Colour | Meaning |
|---|---|---|
| `TRUE` | | the caller reports the variant as written in `FMP` |
| `FALSE` | red | the caller does not report it |
| `DISAGREEMENT` | amber | the caller reports it on the other side of a threshold, for example as major where the other caller has it as minor; `FMP` follows the averaged frequency (`--disagreement_average`) |
| `DISABLED` | grey | Mutect2's minor calls are left out in this region (`--mutect2_disabled_regions`) |
| `LOW`, `OK` | red, none | in a `LOW` row: which caller has too few reads |

The colour of the `FMP` cell shows the kind of row at a glance: red where the
callers do not agree (`FALSE` or `DISAGREEMENT`), yellow for `LOW`, orange for
a point heteroplasmy, blue for a lowercase length variant, and green for
variants both callers report (a lighter green for insertions and deletions).
Red and yellow rows are the ones to look at. The columns that come from
FDSTOOLS have a green background, those from Mutect2 a blue one.

## Main Options

| Option | Default | Meaning |
|---|---|---|
| `--reads` | `raw_data/*_{R1,R2}_001.fastq.gz` | input files |
| `--outdir` | `results` | output folder |
| `--min_vf_FDS`, `--min_vf_MT2` | 5 | a base counts as present from this share (%), per caller |
| `--lh_thresh` | 10 | length variants: under 10% not reported, 10–90% minor (lowercase), above 90% major |
| `--depth` | 10 | amplicons with fewer reads are marked LOW |
| `--frame` | `separate` | how 16180–16193 and 300–315 are written: `separate` or `shared` |
| `--mutect2_disabled_regions` | `none` | complex regions where Mutect2's minor calls are left out: `none`, `all`, or e.g. `300-315,16180-16193` |
| `--disagreement_average` | `plain` | how the two callers' frequencies are averaged when they disagree: `plain` or `depth_weighted` |
| `--mapQ` | 30 | reads with a lower mapping quality are removed |
| `--skip_numt_filter` | `false` | skip rtn (for simulated reads) |
| `--rtn_positions` | `resources/rtn_files/positions/recurrent_positions.tsv` | positions rtn sees at base quality 2 (`none`: off); only an rtn that uses base qualities, such as the fork at <https://github.com/PeterResutik/RtN> (branch `fixes`), counts them less |
| `--publish_bams` | `true` | copy BAM and FASTQ intermediates to the results folder |

All options and their values are printed at the start of each run and saved in
`p00_parameters.txt`.

## Configuration

### Profiles

* `-profile docker`: Runs the pipeline using Docker with all dependencies pre-installed.
* `-profile local`: Runs the pipeline using locally installed tools (requires Conda environment).

## Tests

The Python scripts in `resources/scripts` have unit tests in `tests/`. Run them
from the repository root in the Conda environment:

```bash
conda activate FMP-NimaGen
pytest
```

## Known Limitations

- Mutect2's frequency for a substitution is too high where some reads have no
  base at that position (a substitution and a deletion at one position). In a
  test with 50% rCRS, 30% substituted and 20% deleted, Mutect2 gave the
  deletion 20% but the substitution 42%: it counts the reads without a base
  toward the substitution. Only the frequency shown is affected: whether rCRS
  is still present is decided from Mutect2's read counts, which give rCRS 50%
  there. In a repeat (two equal bases in a row) Mutect2 counts those reads as
  rCRS instead. FDSTOOLS counts each molecule once and is not affected. In the
  C-stretches, where this is most frequent, Mutect2 can be left out with
  `--mutect2_disabled_regions`.

## Cleaning Up

### Remove Cache and Temporary Files

Delete Nextflow's cache and run history (`-resume` no longer works afterwards):

```bash
rm -rf .nextflow
```

Clean up intermediate pipeline files:

```bash
nextflow clean -f
```

### Remove Output Files (Optional)

To delete generated results:

```bash
rm -r work
rm -r results
```

> Use with caution: this permanently deletes output data.

## Citation

If this pipeline is used in research, please cite this repository and the associated tools (e.g., GATK, FDSTools, FLASH, Cutadapt, Nextflow).

## Contributing

Contributions are welcome. Please open an issue or pull request via GitHub if you encounter problems or have suggestions.

## Contact

For support or feedback, submit an issue on the [GitHub repository](https://github.com/PeterResutik/FMP-NimaGen).

## License

FMP-NimaGen is released under the [MIT License](LICENSE).

Third-party material in this repository keeps its own terms:

- `resources/rtn_files/humans/humans_NimaGen_2026-10.fa.bz2` is built from
  [mitoLEAF](https://github.com/forensicgenomics/mitoLeaf) (Mozilla Public
  License 2.0), GenBank, and the 1000 Genomes Project and HGDP data (see
  `resources/rtn_files/humans/build/README.md`).
- `resources/rtn_files/humans/humans.fa.bz2` and `resources/rtn_files/numts/`
  come from [RtN](https://github.com/Ahhgust/RtN).

The tools the pipeline runs (FDSTOOLS, GATK, bwa, samtools, rtn and others) are
not part of this repository and are covered by their own licences.
