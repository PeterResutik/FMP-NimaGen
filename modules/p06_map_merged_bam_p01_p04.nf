process p06_map_merged_bam_p01_p04 {
    tag "p06: bwa mem on $sample_id"
    // publishDir "$params.outdir/p06_mapped_wo_scb_merged_bam/${sample_id}", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged)
    path reference
    path index_files
    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${sample_id}_wo_scb_merged.bam"), path("${sample_id}_wo_scb_merged.bam.bai")

    script:
    """
    bwa mem ${reference} ${fastq_wo_scb_merged} | samtools view -Sb - | samtools sort -o ${sample_id}_wo_scb_sorted.bam

    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${sample_id}_wo_scb_sorted.bam  -o ${sample_id}_wo_scb_merged.bam
    samtools index ${sample_id}_wo_scb_merged.bam
    """
}
