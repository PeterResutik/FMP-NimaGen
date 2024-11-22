nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/data/rCRS.fasta"
// params.multiqc = "$baseDir/multiqc"
params.outdir = "results"

println """\
         m t D N A - N i m a G e n  P I P E L I N E    
         ===================================
         transcriptome: ${params.reference}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

// read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

// read_pairs_ch.view()

process INDEX {   

    input:
    path reference

    output:
    path("${reference}.*")

    script:
    """
    bwa index $reference 
    """
}

process MAPPING {
    tag "bwa mem on $sample_id"

    input:
    path reference
    path index_files
    tuple val(sample_id), path(reads) // Paired-end FASTQ files

    output:
    path "${sample_id}.bam"

    script:
    """
    bwa mem ${index_files[0].baseName} ${reads[0]} ${reads[1]} | samtools view -Sb - > ${sample_id}.bam
    """
}


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    index_ch = INDEX(params.reference)
    mapping_ch = MAPPING(params.reference, index_ch, read_pairs_ch)
}