process p05_trim_merged_fastq_p04 {
    tag "p05: cutadapt on $sample_id"
    // publishDir "$params.outdir/p05_wo_scb_merged_trimmed_fastq/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged)
    path left_primers
    path right_primers

    output:
    tuple val(sample_id), path("${sample_id}_wo_scb_merged_trimmed.fastq")

    script:
    """
    # X instead of ^ and \$ (non-internal): a primer may be incomplete at the read end, as where the first bases of the
    # oligo carry an error that bwa clips and p03 removes, but never inside the read; at least primer_min_overlap of its bases must be there
    tr '^\$' XX < $left_primers > left_primers_noninternal.fasta
    tr '^\$' XX < $right_primers > right_primers_noninternal.fasta
    cutadapt -g file:left_primers_noninternal.fasta -O $params.primer_min_overlap -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_wo_scb_merged_trimmed_left.fastq $fastq_wo_scb_merged
    cutadapt -a file:right_primers_noninternal.fasta -O $params.primer_min_overlap -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_wo_scb_merged_trimmed.fastq ${sample_id}_wo_scb_merged_trimmed_left.fastq
    """
}
