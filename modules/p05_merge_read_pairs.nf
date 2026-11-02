process p05_merge_read_pairs {
    tag "p05: flash on $sample_id"
    // publishDir "$params.outdir/p05_merge_read_pairs/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_r1_wo_scb), path(fastq_r2_wo_scb)
    
    output:
    tuple val(sample_id), path("${sample_id}_wo_scb_merged.extendedFrags.fastq")

    script:
    """
    # -t 1: with more threads FLASH writes the merged reads in a different order on every run
    flash ${fastq_r1_wo_scb} ${fastq_r2_wo_scb} -m $params.min_overlap -M $params.max_overlap -x $params.max_mismatch_density -O -t 1 -o ${sample_id}_wo_scb_merged
    """
}
