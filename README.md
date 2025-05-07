# mtDNA Processing Pipeline

This pipeline processes mitochondrial DNA (mtDNA) sequencing data generated using the Nimagen RC-PCR Kit. It is implemented using Nextflow, enabling streamlined execution across different computational environments.

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

To check your Java version:

```bash
java -version
```

If you don’t have Java 17+, it is recommended to install it via [SDKMAN!](https://sdkman.io/):

#### Install Java via SDKMAN!

1. **Install SDKMAN!**

```bash
curl -s https://get.sdkman.io | bash
```

2. **Restart your terminal**, then run:

```bash
sdk install java 17.0.10-tem
```

3. **Verify installation:**

```bash
java -version
```

#### Install Nextflow

Once Java is set up:

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /usr/local/bin  # Or another directory in your PATH
nextflow -version
```

### Clone the Repository

```bash
git clone https://github.com/PeterResutik/mtDNA-NimaGen.git
cd mtDNA-NimaGen
```

### Organize Input Data

Move raw sequencing files into the designated folder:

```bash
cp /path/to/your/FASTQ/*fastq.gz raw_data/
```

### Create Conda Environment

> **Note**: This pipeline assumes you have [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed.
> If not, download and install Miniconda from the [official website](https://docs.conda.io/en/latest/miniconda.html) before proceeding.

Create the Conda environment from the provided YAML file:

### Prepare Reference Files

Navigate to the reference files directory and index the genome:

```bash
cd resources/rtn_files
bunzip2 humans.fa.bz2
bwa index humans.fa
```
> ** Note for Docker users**  
> These commands must be run **outside** the pipeline using locally installed tools.  
> The Docker profile only affects tools used *during* pipeline execution.  
> You must have `bwa` and `bunzip2` installed on your system to prepare the reference genome.
>
>   You can install `bwa` easily using Conda:  
> ```bash
> conda install -c bioconda bwa
> ```

## Running the Pipeline

### Execution Options

The pipeline can be executed using either Docker or local dependencies.

#### Running with Docker (Recommended)

```bash
nextflow run main.nf -profile docker
```

This method ensures that all required dependencies are available.

#### Resuming a Previous Run

To avoid re-processing completed steps, use:

```bash
nextflow run main.nf -profile docker -resume
```

## Local Installation

To run the pipeline locally, create the Conda environment from the provided YAML file:

```bash
conda env create -f mtDNA_env.yml
```

<!-- Additional Python dependencies:

```bash
conda install -n mtDNA_env -c conda-forge pandas biopython matplotlib
``` -->

### Virtual Environment

A Conda environment is recommended to isolate dependencies. Before running the pipeline, activate it:

```bash
conda activate mtDNA_env

nextflow run main.nf -profile local
```

## Configuration

### Profiles

- `-profile docker`: Uses a Docker container with pre-installed dependencies.
- `-profile local`: Uses software installed locally.

## Cleaning Up

### Remove Cache and Temporary Files

To start from a clean state, delete cached data:

```bash
rm -rf .nextflow.cache
```

To remove intermediate work files:

```bash
nextflow clean -f
```

### Remove Output Files (Optional)

To delete previous results:

```bash
rm -r work
rm -r results
```

Use this with caution, as it permanently deletes output data.

## Citation

If this pipeline is used in research, please cite this repository along with the associated tools.

## Contributing

Contributions are welcome. To report issues or suggest improvements, open an issue or pull request on GitHub.

## Contact

For support, submit an issue on the GitHub repository.


