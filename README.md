# mtDNA Processing Pipeline

Welcome to the **mtDNA Processing Pipeline**, an optimized and reproducible workflow for processing raw FASTQ files from mitochondrial DNA (mtDNA) sequences generated using the **Nimagen RC-PCR Kit**. This pipeline is implemented using [Nextflow](https://www.nextflow.io/), ensuring seamless execution across various computational environments.

---

## Overview

The pipeline automates essential steps in mtDNA sequencing data analysis, including:

1. **Soft-clip Removal**: Eliminates soft-clipped bases for improved read quality.
2. **Read Merging**: Uses FLASH to merge overlapping paired-end reads for enhanced sequence accuracy.
3. **Primer Trimming**: Utilizes Cutadapt to remove primers and unwanted sequences.
4. **Variant Calling**: Employs GATK's Mutect2 to detect mtDNA variants.
5. **Forensic Nomenclature Conversion**: Converts variants into a forensic nomenclature consistent with the **EMPOP** database.

This pipeline provides a **transparent**, **reproducible**, and **scalable** workflow for scientists performing mtDNA analyses.

<!-- ---

## Features

- **Reproducible**: Ensures consistent results across multiple runs.
- **Scalable**: Supports execution on local machines, clusters, and cloud environments.
- **Customizable**: Allows fine-tuning of parameters for read merging, trimming, and variant calling. -->

---

## Requirements

### Software

- **Nextflow**: >= 22.10.0
- **Java**: >= 8
- **Docker** or **Singularity** (optional, for containerized execution)

### Dependencies

The pipeline relies on the following tools, which can be installed via conda or containerized environments:

- [FLASH (for read merging)](https://ccb.jhu.edu/software/FLASH/)
- [Cutadapt (for primer trimming)](https://cutadapt.readthedocs.io/)
- [BWA (for sequence alignment)](http://bio-bwa.sourceforge.net/)
- [GATK Mutect2 (for variant calling)](https://gatk.broadinstitute.org/)

---

## Installation

Clone the repository and navigate to the pipeline directory:

```bash
git clone https://github.com/<your-username>/mtDNA-pipeline.git
cd mtDNA-pipeline
```

---

## Execution

To run the pipeline, execute the following command:

```bash
nextflow run main.nf
```

---

## Configurable Parameters

The following parameters can be adjusted in the pipeline:

### Read Merging (FLASH)

```bash
params.min_overlap = 10  # Default: 10
params.max_overlap = 140 # Default in FLASH: 65
params.max_mismatch_density = 0.25 # Default in FLASH: 0.25
```

### Primer Trimming (Cutadapt)

```bash
params.quality_cutoff = 25
params.minimum_length = 60
params.maximum_length = 300
```

### Variant Calling (GATK Mutect2)

```bash
params.detection_limit = 0.075
params.mapQ = 30
params.baseQ = 32
params.alignQ = 30
```

---

## Output

- **Processed Reads**: Trimmed and merged reads aligned to the reference genome.
- **Variant Calls**: Identified variants in VCF format.
- **Forensic Format Variants**: Variants converted to forensic nomenclature compatible with EMPOP.

---

## Citation

If you use this pipeline in your research, please cite this repository and the tools used.

---

For any issues or questions, feel free to open an issue on the GitHub repository.

Happy Analyzing!

