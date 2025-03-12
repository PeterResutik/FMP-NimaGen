# mtDNA Processing Pipeline

This pipeline processes mitochondrial DNA (mtDNA) sequencing data generated using the Nimagen RC-PCR Kit. It is implemented using Nextflow, enabling streamlined execution across different computational environments.

## Overview

The pipeline automates the following steps:

1. **Soft-clip Removal**: Removes soft-clipped bases for improved read quality.
2. **Read Merging**: Uses FLASH to merge overlapping paired-end reads.
3. **Primer Trimming**: Uses Cutadapt to remove primers and unwanted sequences.
4. **Variant Calling**: Uses GATK Mutect2 to detect mtDNA variants.
5. **Forensic Nomenclature Conversion**: Converts variants into a nomenclature consistent with the EMPOP database.

The workflow is designed for reproducibility and scalability.

## Installation and Setup

### Clone the Repository

```bash
git clone https://github.com/PeterResutik/mtDNA-NimaGen.git
cd mtDNA-NimaGen
```

### Prepare Reference Files

Navigate to the reference files directory and index the genome:

```bash
cd resources/rtn_files
bunzip2 humans.fa.bz2 && bwa index humans.fa
```

### Organize Input Data

Move raw sequencing files into the designated folder:

```bash
mkdir -p raw_data
cp /path/to/your/FASTQ/*fastq.gz raw_data/
```

## Running the Pipeline

### Execution Options

The pipeline can be executed using either Docker or local dependencies.

#### Running with Docker (Recommended)

```bash
nextflow run main.nf -profile docker
```

This method ensures that all required dependencies are available.

#### Running Locally

If running locally, first install dependencies as outlined below. Then execute:

```bash
nextflow run main.nf -profile local
```

#### Resuming a Previous Run

To avoid re-processing completed steps, use:

```bash
nextflow run main.nf -profile docker -resume
```

## Local Installation

To run the pipeline locally, install the required software using Conda:

```bash
conda create -n mtDNA_env -c bioconda -c conda-forge -y \
    nextflow=24.10.0 \
    python=3.9 \
    cutadapt=4.6 \
    samtools=1.21 \
    bwa=0.7.18 \
    flash=1.2.11 \
    csvtk=0.31.0 \
    fastqc=0.11.9 \
    bcftools=1.11-35-g8a744dd \
    gatk4=4.6.1.0
```

Additional Python dependencies:

```bash
conda install -n mtDNA_env -c conda-forge pandas biopython matplotlib
```

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

### Configurable Parameters

Adjust parameters in `nextflow.config` or specify them at runtime.

#### Read Merging (FLASH)

```bash
params.min_overlap = 10  
params.max_overlap = 140 
params.max_mismatch_density = 0.25  
```

#### Primer Trimming (Cutadapt)

```bash
params.quality_cutoff = 25
params.minimum_length = 60
params.maximum_length = 300
```

#### Variant Calling (GATK Mutect2)

```bash
params.detection_limit = 0.075
params.mapQ = 30
params.baseQ = 32
params.alignQ = 30
```

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
rm -r work.nosync
rm -r results.nosync
```

Use this with caution, as it permanently deletes output data.

## Output Files

The pipeline generates the following:

| File Type | Description |
|-----------|-------------|
| Processed Reads | Trimmed and merged sequences aligned to the reference genome |
| Variant Calls | Identified variants in VCF format |
| Forensic Format Variants | Variants formatted for EMPOP |

## Citation

If this pipeline is used in research, please cite this repository along with the associated tools.

## Contributing

Contributions are welcome. To report issues or suggest improvements, open an issue or pull request on GitHub.

## Contact

For support, submit an issue on the GitHub repository.
