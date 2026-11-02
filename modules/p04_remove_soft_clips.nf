process p04_remove_soft_clips {
    tag "p04: removing scb from $sample_id"
    publishDir { "$params.outdir/p04_remove_soft_clips/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*sorted.bam*', enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(sam_r1), path(sam_r2)
    path python_script_remove_scb
    path reference
    path index_files
    path left_primers
    path right_primers

    output:
    tuple val(sample_id), path("${sample_id}_R1_wo_scb.bam"), path("${sample_id}_R2_wo_scb.bam")
    tuple val(sample_id), path("${sample_id}_R1_wo_scb_sorted.bam"), path("${sample_id}_R2_wo_scb_sorted.bam"), path("${sample_id}_R1_wo_scb_sorted.bam.bai"), path("${sample_id}_R2_wo_scb_sorted.bam.bai")
    tuple val(sample_id), path("${sample_id}_R1_R2_wo_scb_sorted.bam"), path("${sample_id}_R1_R2_wo_scb_sorted.bam.bai")
    tuple val(sample_id), path("${sample_id}_R1_wo_scb.fastq"), path("${sample_id}_R2_wo_scb.fastq"), emit: fastq_wo_scb_ch


    script:
    """ 
    # soft-clipped bases go, and so do the few bases a read reaches past a primer boundary (read-through)
    python $python_script_remove_scb ${sam_r1} ${sample_id}_R1_wo_scb.sam ${reference} ${left_primers} ${right_primers}
    python $python_script_remove_scb ${sam_r2} ${sample_id}_R2_wo_scb.sam ${reference} ${left_primers} ${right_primers}

    samtools view -Sb ${sample_id}_R1_wo_scb.sam > ${sample_id}_R1_wo_scb.bam
    samtools view -Sb ${sample_id}_R2_wo_scb.sam > ${sample_id}_R2_wo_scb.bam

    samtools fastq ${sample_id}_R1_wo_scb.bam > ${sample_id}_R1_wo_scb.fastq    
    samtools fastq ${sample_id}_R2_wo_scb.bam > ${sample_id}_R2_wo_scb.fastq   

    bwa mem ${reference} ${sample_id}_R1_wo_scb.fastq ${sample_id}_R2_wo_scb.fastq  > ${sample_id}_R1_R2_wo_scb.sam  
    samtools view -bS ${sample_id}_R1_R2_wo_scb.sam | samtools sort -o ${sample_id}_R1_R2_wo_scb_sorted.bam
    samtools index ${sample_id}_R1_R2_wo_scb_sorted.bam

    samtools sort -o ${sample_id}_R1_wo_scb_sorted.bam ${sample_id}_R1_wo_scb.bam
    samtools sort -o ${sample_id}_R2_wo_scb_sorted.bam ${sample_id}_R2_wo_scb.bam
    samtools index ${sample_id}_R1_wo_scb_sorted.bam
    samtools index ${sample_id}_R2_wo_scb_sorted.bam
    """
}
