# Forensic mtDNA Pipeline (FMP)

This pipeline processes mitochondrial DNA (mtDNA) sequencing data generated using the Nimagen RC-PCR Kit. It is implemented using Nextflow, enabling streamlined execution across different computational environments.

## Table of Contents

- [Overview](#overview)  
- [Installation and Setup](#installation-and-setup)  
  - [Nextflow Requirements and Installation](#nextflow-requirements-and-installation)  
  - [Clone the Repository](#clone-the-repository)  
  - [Organize Input Data](#organize-input-data)  
- [Running the Pipeline](#running-the-pipeline)  
  - [Running with Docker (Recommended)](#running-with-docker-recommended)  
  - [Running Locally with Conda](#running-locally-with-conda)  
- [Configuration](#configuration)  
- [Cleaning Up](#cleaning-up)  
- [Citation](#citation)  
- [Contributing](#contributing)  
- [Contact](#contact)  

## Overview

The pipeline automates the following steps:

1. **Soft-clip Removal**: Removes soft-clipped bases to improve alignment quality.  
2. **Read Merging**: Uses FLASH to merge overlapping paired-end reads.  
3. **Primer Trimming**: Uses Cutadapt to remove primers and unwanted flanking sequences.  
4. **Variant Calling**: Uses both **GATK Mutect2** and **FDSTools** to call mtDNA variants.  
5. **Variant Comparison**: Merges and compares outputs from both tools to support interpretation and distinguish true variants from false positives.  

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
git clone https://github.com/PeterResutik/mtDNA-NimaGen.git
cd mtDNA-NimaGen
```

### Organize Input Data

Place your raw sequencing files into the expected directory:

```bash
cp /path/to/your/FASTQ/*fastq.gz raw_data/
```

## Running the Pipeline

### Running with Docker (Recommended)

This is the recommended and most portable way to run the pipeline:

```bash
nextflow run main.nf -profile docker
```

To resume a previous run and skip already completed steps:

```bash
nextflow run main.nf -profile docker -resume
```

> **Note**: You do **not** need to install Conda or any dependencies if using Docker.

### Running Locally with Conda

If you prefer to run the pipeline with your own local environment:

1. Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Conda](https://docs.conda.io/en/latest/) installed.
2. Create the environment:

```bash
conda env create -f FMP-NimaGen.yml
```

3. Activate the environment and run the pipeline:

```bash
conda activate FMP-NimaGen.yml
nextflow run main.nf -profile local
```

## Configuration

### Profiles

* `-profile docker`: Runs the pipeline using Docker with all dependencies pre-installed.
* `-profile local`: Runs the pipeline using locally installed tools (requires Conda environment).

## Cleaning Up

### Remove Cache and Temporary Files

Delete Nextflow's execution cache:

```bash
rm -rf .nextflow.cache
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

For support or feedback, submit an issue on the [GitHub repository](https://github.com/PeterResutik/mtDNA-NimaGen).

