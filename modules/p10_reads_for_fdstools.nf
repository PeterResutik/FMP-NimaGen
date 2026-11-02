process p10_reads_for_fdstools {
    tag "p10: reads rtn kept on $sample_id"
    publishDir { "$params.outdir/p10_reads_for_fdstools/${sample_id}" }, mode: 'copy', overwrite: true, enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(bam_wo_scb_merged), path(bam_wo_scb_merged_index), path(kept_names)

    output:
    tuple val(sample_id), path("${bam_wo_scb_merged.baseName}.rtn.bam"), path("${bam_wo_scb_merged.baseName}.rtn.bam.bai"), path("${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq")

    script:
    // the NUMT decision comes from p09, where rtn judges the same reads without their primers: the primers
    // are the kit's oligos, and where they carry a rare rCRS base (750A, 8860A) every read would differ from almost every genome
    """
    samtools view -b -h -q $params.mapQ -N ${kept_names} ${bam_wo_scb_merged} > ${bam_wo_scb_merged.baseName}.rtn.bam
    samtools index ${bam_wo_scb_merged.baseName}.rtn.bam
    samtools fastq ${bam_wo_scb_merged.baseName}.rtn.bam > ${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq
    """
}
