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
5. **NUMT filtering** (`p08`–`p09`): rtn removes reads that no human mitogenome in `humans_NimaGen.fa` explains closely enough, which removes reads from nuclear copies of mtDNA (NUMTs); reads below `--mapQ` are removed too.  
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

- `resources/rtn_files/humans/humans_NimaGen.fa.bz2` contains sequences derived
  from [mitoLEAF](https://github.com/forensicgenomics/mitoLeaf), which is licensed
  under the Mozilla Public License 2.0 (see
  `resources/rtn_files/humans/build/README.md`).
- `resources/rtn_files/humans/humans.fa.bz2` and `resources/rtn_files/numts/`
  come from [RtN](https://github.com/Ahhgust/RtN).

The tools the pipeline runs (FDSTOOLS, GATK, bwa, samtools, rtn and others) are
not part of this repository and are covered by their own licences.
