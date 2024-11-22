#!/bin/bash -ue
bwa mem rCRS.fasta Ctr-1-10_S11_L001_R1_001.fastq.gz Ctr-1-10_S11_L001_R2_001.fastq.gz | samtools view -Sb - > Ctr-1-10_S11_L001.bam
