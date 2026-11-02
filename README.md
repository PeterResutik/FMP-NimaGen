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
  - [Running with Docker or Singularity](#running-with-docker-or-singularity)  
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

1. **Mapping and soft-clip removal** (`p03`–`p04`): read pairs are mapped to rCRS with bwa, and soft-clipped bases are removed.  
2. **Read merging** (`p05`): FLASH merges each read pair into one read.  
3. **Primer trimming** (`p06`): Cutadapt removes the primers and drops reads without a primer or outside the length limits. The primer files are described in `resources/primers/README.md`.  
4. **Mapping** (`p07`–`p08`): the merged reads are mapped to rCRS with their primers (for FDSTOOLS) and without (for rtn).  
5. **NUMT filtering** (`p09`–`p10`): rtn removes reads that no human mitogenome in `humans_NimaGen_2026-10-30.fa` explains closely enough (`--rtn_likelihood`), which removes reads from nuclear copies of mtDNA (NUMTs); reads below `--mapQ` are removed too. rtn judges the reads without their primers (`p09`); `p10` keeps the same reads, by name, for FDSTOOLS. At the 50 positions that change most often in human mtDNA (`resources/rtn_files/positions/`), rtn sees the reads at base quality 2, so that differences there count little; FDSTOOLS keeps the original qualities.  
6. **Quality control** (`p11`): read depth per amplicon.  
7. **Variant calling** (`p12`): **FDSTOOLS** on the reads with primers. It distrusts bases below `--baseQ` (30): a read that differs from a more common sequence of the same length only at such bases counts for that sequence.  
8. **Report** (`p13`): FDSTOOLS' calls are written in one notation into one Excel table per sample.  

Each step is a module in `modules/` (one Nextflow process per file), numbered in the order the steps run (`p00` to `p13`); `main.nf` connects them, and each result folder has the name of its step.

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

### Install rtn

The NUMT filter runs [rtn](https://github.com/Ahhgust/RtN), which is not part of the
conda environment. Build it from the fork with three fixes to how rtn scores reads
(it uses base qualities, which `--rtn_positions` needs, and no longer drops reads
with length variants in the C-stretches); this is tested on Linux:

```bash
git clone -b fixes https://github.com/PeterResutik/RtN.git
cd RtN
conda create -n rtn-build -c conda-forge jsoncpp
conda activate rtn-build
PREFIX=$CONDA_PREFIX ./build.sh
cp rtn /usr/local/bin/   # or any directory in your PATH
```

`build.sh` also needs git, make, autoconf, automake, a C/C++ compiler and the
development files of zlib, bzip2 and xz. It does not yet build on Apple Silicon.
rtn as released also runs with the pipeline, but without the fixes: it ignores base
qualities, so `--rtn_positions` has no effect, and it drops some reads with length
variants.

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
Nextflow runs as many tasks at once as the computer has cores; to keep the
computer responsive during a run, limit them with `-qs`, e.g. `-qs 4`.

### Running with Docker or Singularity

The image `peterresutik/nimagen-pipeline:0.5.0` holds every tool the pipeline
runs (the Conda environment and rtn with its three fixes); Nextflow, Java and the
pipeline folder stay outside it. FDSTOOLS downloads the chrM reference sequence
from Ensembl the first time it runs and keeps it in `~/.strnaming-cache`; the
image already holds it, so the containers need neither internet access nor a
home folder. With Docker:

```bash
nextflow run main.nf -profile docker
```

With Singularity (or Apptainer), for example on a cluster without Docker:

```bash
nextflow run main.nf -profile singularity
```

Nextflow pulls the image on the first run (Singularity converts it and keeps it in
`work/singularity/`, or in `NXF_SINGULARITY_CACHEDIR` if set). The image is built
for x86-64 Linux, because rtn's alignment library uses x86 instructions; on Apple
Silicon, Docker runs it through emulation, which is slower but gives the same
results. To build the image yourself:

```bash
docker build --platform linux/amd64 -t peterresutik/nimagen-pipeline:0.5.0 .
```

## Input and Output

Read files must be named `<sample>_R1_001.fastq.gz` and `<sample>_R2_001.fastq.gz`
and are taken from `raw_data/` by default. Results go to `results/`:

| Folder | Contents |
|---|---|
| `p00_parameters.txt` | the options of the run |
| `p11_read_depth/<sample>/` | read depth per amplicon (plot and tables) |
| `p12_fdstools/<sample>/` | FDSTOOLS files (sequences per amplicon, statistics, HTML) |
| `p13_report/<sample>_separate_frame.xlsx` | **the report** (`_shared_frame` with `--frame shared`) |

BAM and FASTQ intermediates are also copied unless `--publish_bams false` is set.

## Reading the Report

The report has one row per variant, sorted by position. The first column,
`FMP`, is the call. The columns after it show the reads FDSTOOLS counted for it.

### Notation

| Example | Meaning |
|---|---|
| `A263G` | substitution: rCRS base, position, base found |
| `A523-` | deletion of the rCRS base at 523 |
| `-309.1C` | insertion of a C after 309; a second inserted base is `-309.2C` |
| `T195Y` | point heteroplasmy: two or more bases at 195, written with the IUPAC code (`Y` is C and T) |
| `-309.2c` | lowercase: a length variant on part of the molecules only (between `--lh_thresh` and 100 minus `--lh_thresh` percent) |
| `LOW` | not a variant: an amplicon with fewer reads than `--depth`, named in `marker`; calls in `marker_range` rest on few reads |

In 16180–16193 and 300–315 the labels follow `--frame`: `separate` writes
each run of equal bases by its length, `shared` lays the sequence from the
left, as mitoLEAF and EMPOP do. In both frames, one base inside an A-run that
has all its A's is an insertion: AAACA at 16180 is `-16182.1C`, where mitoLEAF
writes 16184A. Outside these regions a molecule is written with the fewest
changes from rCRS. That can differ from the spelling of the phylogeny: L1c1c
and L5b1 are written `G247-` where the phylogeny writes `247A 249del`. The
molecule is the same.

### Columns

| Column | Meaning |
|---|---|
| `FMP` | the call |
| `vf_FDS` | share of the reads with the variant, in percent; with several bases at one position, one value per base (`C 40, A 10`) |
| `rd_FDS`, `rd_total_FDS` | reads with the variant, and all reads at the position; in a `LOW` row, the reads of the amplicon |
| `variant_note` | remarks; not written with `--variant_note false` |
| `dominant_molecule` | in 16180–16193, 300–315 and 57–60: what the molecule with the most reads carries at this position, with that molecule's share of the reads, e.g. `-16193.1C (40.9%)`; empty if it carries no change there |
| `marker`, `marker_range`, `num_markers` | the amplicons FDSTOOLS' depth and frequency come from (every amplicon with reads covering the position), their ranges, and their number |

FDSTOOLS gives the reads and the share of the variant (`rd_FDS`, `vf_FDS`);
the reads without it follow from `rd_total_FDS`. Reads that run from one
amplicon's left primer to the right primer of the next amplicon in the same
primer pool (001 to 003, 002 to 004, 098 to 100) are whole mtDNA molecules,
but FDSTOOLS cannot assign them to one amplicon, so they do not count towards
its depth.

### Colours

The colour of the `FMP` cell shows the kind of row at a glance: yellow for
`LOW`, orange for a point heteroplasmy, blue for a lowercase length variant,
and green for the other variants (a lighter green for insertions and
deletions). Yellow rows are the ones to look at. The other columns have a
green background.

## Main Options

| Option | Default | Meaning |
|---|---|---|
| `--reads` | `raw_data/*_{R1,R2}_001.fastq.gz` | input files |
| `--outdir` | `results` | output folder |
| `--min_vf_FDS` | 2 | a base counts as present from this share (%) |
| `--max_overlap` | 226 | FLASH counts the mismatches of an overlap over at most this many bases; set it to the read length, so that the whole overlap is checked |
| `--lh_thresh` | 10 | length variants: under 10% not reported, 10–90% minor (lowercase), above 90% major |
| `--depth` | 20 | amplicons with fewer reads are marked LOW |
| `--frame` | `separate` | how 16180–16193 and 300–315 are written: `separate` or `shared` |
| `--variant_note` | `true` | `false` leaves the `variant_note` column out of the report |
| `--mapQ` | 30 | reads with a lower mapping quality are removed |
| `--skip_numt_filter` | `false` | skip rtn (for simulated reads) |
| `--rtn_positions` | `resources/rtn_files/positions/recurrent_positions.tsv` | positions rtn sees at base quality 2 (`none`: off); only an rtn that uses base qualities, such as the fork at <https://github.com/PeterResutik/RtN> (branch `fixes`), counts them less |
| `--rtn_likelihood` | `1e-7` | the minimum likelihood a read must have under the human mitogenomes (rtn `-L`; rtn's own default is 1e-5); a lower value removes fewer reads |
| `--rtn_memory` | `4 GB` | memory per rtn task (`p09`), which holds the humans index (about 2 GB); Nextflow starts no more tasks than the computer's memory holds |
| `--num_threads` | 6 | threads of FDSTOOLS `tssv` (`p12`); Nextflow reserves as many cores for each task, at most as many as the computer has |
| `--publish_bams` | `true` | copy BAM and FASTQ intermediates to the results folder |

All options and their defaults are in the `params` block of `nextflow.config`. A
value given on the command line (`--name value`) or in a file passed with
`-params-file` overrides the default. All options and their values are printed at
the start of each run and saved in `p00_parameters.txt`.

## Configuration

### Profiles

* `-profile local`: runs the tools installed in the active Conda environment (`FMP-NimaGen.yml`) and rtn.
* `-profile docker`: runs every tool in the Docker image (see Running with Docker or Singularity).
* `-profile singularity`: the same image with Singularity or Apptainer.

## Tests

The Python scripts in `resources/scripts` have unit tests in `tests/`. Run them
from the repository root in the Conda environment:

```bash
conda activate FMP-NimaGen
pytest
```

## Known Limitations

- A sequencing error that comes back at the same position in many samples can
  reach `--min_vf_FDS` and appear as a heteroplasmy. Before FDSTOOLS distrusted
  bases below `--baseQ` (0.3.0), this happened in the samples tested at 3407,
  3413 and 3460 (amplicon 023, whose R2 loses quality after its first cycles),
  at 8541 (amplicon 053) and at 10231 and 10266 (amplicon 063), at up to 4%.
  A higher detection threshold (`--min_vf_FDS`) leaves them out.

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

If this pipeline is used in research, please cite this repository and the associated tools (e.g., FDSTools, FLASH, Cutadapt, Nextflow).

## Contributing

Contributions are welcome. Please open an issue or pull request via GitHub if you encounter problems or have suggestions.

## Contact

For support or feedback, submit an issue on the [GitHub repository](https://github.com/PeterResutik/FMP-NimaGen).

## License

FMP-NimaGen is released under the [MIT License](LICENSE).

Third-party material in this repository keeps its own terms:

- `resources/rtn_files/humans/humans_NimaGen_2026-10-30.fa.bz2` (and the earlier
  `humans_NimaGen_2026-10-19.fa.bz2`) are built from
  [mitoLEAF](https://github.com/forensicgenomics/mitoLeaf) (Mozilla Public
  License 2.0), GenBank, and the 1000 Genomes Project and HGDP data (see
  `resources/rtn_files/humans/build/README.md`).
- `resources/rtn_files/numts/` comes from [RtN](https://github.com/Ahhgust/RtN).

The tools the pipeline runs (FDSTOOLS, bwa, samtools, rtn and others) are
not part of this repository and are covered by their own licences.
