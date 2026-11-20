process p09_filter_numts {
    tag "p09: rtn on $sample_id"
    memory params.rtn_memory
    publishDir { "$params.outdir/p09_filter_numts/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*.bam*', enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(bam_wo_scb_merged_trimmed), path(bam_wo_scb_merged_trimmed_index), path(read_depth_txt)
    path amplicon_middle_positions
    path humans_index
    val humans_base
    path numts_index
    val numts_base
    val rtn_lowq
    path python_script_rtn_qualities
    path python_script_mark_removed_reads


    output:
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam"), path("${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam.bai"), path(read_depth_txt), path("${bam_wo_scb_merged_trimmed.baseName}_read_depth_wo_NUMTs.txt"), emit: bams
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}_kept_names.txt"), path("${bam_wo_scb_merged_trimmed.baseName}_removed_reads.tsv"), emit: read_names
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}_all_reads.rtn.bam"), path("${bam_wo_scb_merged_trimmed.baseName}_all_reads.rtn.bam.bai"), emit: all_reads

    script:
    // rtn -i: indels do not count toward a read's distance, as when rtn picks the closest genome
    // rtn -L: the minimum likelihood a read must have under the human mitogenomes (--rtn_likelihood)
    def base = bam_wo_scb_merged_trimmed.baseName
    def rtn = "rtn -i -L ${params.rtn_likelihood} -h \"${humans_index}/${humans_base}\" -n \"${numts_index}/${numts_base}\" -b"
    // rtn_lowq: rtn evaluates a copy of the reads with base quality 2 at these positions; its decision
    // (MAPQ, tags) is kept and the reads get their own qualities back
    def numt_filter = params.skip_numt_filter.toString().toBoolean() ?
        "cp ${bam_wo_scb_merged_trimmed} ${base}.rtn.bam && samtools index ${base}.rtn.bam" :
        rtn_lowq == 'none' ? "${rtn} ${bam_wo_scb_merged_trimmed}" :
        """set -o pipefail
    samtools view -h ${bam_wo_scb_merged_trimmed} | python ${python_script_rtn_qualities} lower ${rtn_lowq} > ${base}.lowq.sam
    samtools view -b -o ${base}.lowq.bam ${base}.lowq.sam && samtools index ${base}.lowq.bam
    ${rtn} ${base}.lowq.bam
    samtools view -h ${base}.lowq.rtn.bam | python ${python_script_rtn_qualities} restore ${bam_wo_scb_merged_trimmed} > ${base}.restored.sam
    samtools view -b -o ${base}.rtn.bam ${base}.restored.sam && samtools index ${base}.rtn.bam
    rm -f ${base}.lowq.* ${base}.restored.sam"""
    // reads that the mapping (p08) could not align end to end go too: rtn judges only the aligned part of a read (it
    // cuts soft-clipped bases off first); a split read is clipped in every record, so all its parts go
    def clipped = params.skip_numt_filter.toString().toBoolean() ? "" : "| awk -F'\\t' '/^@/ || \$6 !~ /[SH]/' "
    def clip_rule = params.skip_numt_filter.toString().toBoolean() ? 0 : 1
    // all_reads: every read of this step, the removed ones at MAPQ 0, each with its reason in the tag ZF
    """
    ${numt_filter}

    samtools view -h -q $params.mapQ ${bam_wo_scb_merged_trimmed.baseName}.rtn.bam ${clipped}| samtools view -b -o ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam
    samtools index ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam
    samtools view ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam | cut -f1 | sort -u > ${bam_wo_scb_merged_trimmed.baseName}_kept_names.txt
    samtools depth -a -b $amplicon_middle_positions ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam > ${bam_wo_scb_merged_trimmed.baseName}_read_depth_wo_NUMTs.txt

    set -o pipefail
    samtools view -h ${bam_wo_scb_merged_trimmed} | python ${python_script_mark_removed_reads} p09 ${base}.rtn.bam $params.mapQ ${clip_rule} ${base}_removed_reads.tsv \\
        | samtools view -b -o ${base}_all_reads.rtn.bam
    samtools index ${base}_all_reads.rtn.bam
    [ \$(samtools view -c -d ZF:kept ${base}_all_reads.rtn.bam) -eq \$(samtools view -c ${base}_filtered.rtn.bam) ] || { echo "p09: the reads marked kept are not the reads kept" >&2; exit 1; }
    """
}
