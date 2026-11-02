process p03_map_read_pairs {
    tag "p03: bwa mem on $sample_id"
    publishDir { "$params.outdir/p03_map_read_pairs/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*.bam*', enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(reads)
    path reference
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_R1.sam"), path("${sample_id}_R2.sam"), emit: raw_sam_ch
    tuple path("${sample_id}_R1.bam"), path("${sample_id}_R2.bam"), path("${sample_id}_R1.bam.bai"), path("${sample_id}_R2.bam.bai")
    tuple path("${sample_id}_R1_R2.bam"), path("${sample_id}_R1_R2.bam.bai")

    script:
    """
    mv ${reads[0]} tmp.fastq.gz
    # -O: only the whole sequence counts; by default 3 matching bases at the end of any read are cut off
    # --action=retain: the anchor stays in R1, only the bases after it go
    cutadapt -a ${params.adapter} -O ${params.adapter.size()} --action=retain -o ${reads[0]}  tmp.fastq.gz 

    bwa mem ${reference} ${reads[0]} > ${sample_id}_R1.sam 
    bwa mem ${reference} ${reads[1]} > ${sample_id}_R2.sam 
    bwa mem ${reference} ${reads[0]} ${reads[1]} > ${sample_id}_R1_R2.sam  
    
    samtools view -bS ${sample_id}_R1_R2.sam | samtools sort -o ${sample_id}_R1_R2.bam
    samtools index ${sample_id}_R1_R2.bam

    samtools view -bS ${sample_id}_R1.sam | samtools sort -o ${sample_id}_R1.bam
    samtools view -bS ${sample_id}_R2.sam | samtools sort -o ${sample_id}_R2.bam
    samtools index ${sample_id}_R1.bam
    samtools index ${sample_id}_R2.bam 
    """
}
