process p10_reads_for_fdstools {
    tag "p10: reads rtn kept on $sample_id"
    publishDir { "$params.outdir/p10_reads_for_fdstools/${sample_id}" }, mode: 'copy', overwrite: true, enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(bam_wo_scb_merged), path(bam_wo_scb_merged_index), path(kept_names), path(removed_reads)
    path python_script_mark_removed_reads

    output:
    tuple val(sample_id), path("${bam_wo_scb_merged.baseName}.rtn.bam"), path("${bam_wo_scb_merged.baseName}.rtn.bam.bai"), path("${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq"), emit: reads
    tuple val(sample_id), path("${bam_wo_scb_merged.baseName}_all_reads.rtn.bam"), path("${bam_wo_scb_merged.baseName}_all_reads.rtn.bam.bai"), emit: all_reads

    script:
    // the NUMT decision comes from p09, where rtn judges the same reads without their primers: the primers
    // are the kit's oligos, and where they carry a rare rCRS base (750A, 8860A) every read would differ from almost every genome
    // all_reads: every read with its primers, the removed ones at MAPQ 0, each with its reason (from p09) in the tag ZF
    def base = bam_wo_scb_merged.baseName
    """
    samtools view -b -h -q $params.mapQ -N ${kept_names} ${bam_wo_scb_merged} > ${bam_wo_scb_merged.baseName}.rtn.bam
    samtools index ${bam_wo_scb_merged.baseName}.rtn.bam
    samtools fastq ${bam_wo_scb_merged.baseName}.rtn.bam > ${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq

    set -o pipefail
    samtools view -h ${bam_wo_scb_merged} | python ${python_script_mark_removed_reads} p10 ${kept_names} ${removed_reads} $params.mapQ \\
        | samtools view -b -o ${base}_all_reads.rtn.bam
    samtools index ${base}_all_reads.rtn.bam
    [ \$(samtools view -c -d ZF:kept ${base}_all_reads.rtn.bam) -eq \$(samtools view -c ${base}.rtn.bam) ] || { echo "p10: the reads marked kept are not the reads kept" >&2; exit 1; }
    """
}
