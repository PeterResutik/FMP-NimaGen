process p01_index_reference_fasta {
    tag "p01: bwa index on $reference"   
    publishDir "$params.outdir/p01_index", mode: 'copy', overwrite: true

    input:
    path reference

    output:
    path("${reference}.*"), emit: reference

    script:
    """
    bwa index $reference 
    
    samtools faidx $reference 
    """
}
