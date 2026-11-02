process p07_map_merged_trimmed_bam_p01_p05 {
    tag "p07: bwa mem on $sample_id"
    // publishDir "$params.outdir/p07_mapped_wo_scb_merged_trimmed_bam/${sample_id}", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged_trimmed)
    path reference
    path index_files

    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${fastq_wo_scb_merged_trimmed.baseName}.bam"), path("${fastq_wo_scb_merged_trimmed.baseName}.bam.bai"), path("${fastq_wo_scb_merged_trimmed.baseName}_read_depth.txt")

    script:
    """
    # read ends are amplicon ends after primer trimming: clip only clearly wrong ends
    bwa mem -L ${params.clipping_penalty} ${reference} ${fastq_wo_scb_merged_trimmed} | samtools view -Sb - | samtools sort -o ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam

    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam  -o ${fastq_wo_scb_merged_trimmed.baseName}.bam
    samtools index ${fastq_wo_scb_merged_trimmed.baseName}.bam

    samtools depth -a -b $amplicon_middle_positions ${fastq_wo_scb_merged_trimmed.baseName}.bam > ${fastq_wo_scb_merged_trimmed.baseName}_read_depth.txt
    """
}
