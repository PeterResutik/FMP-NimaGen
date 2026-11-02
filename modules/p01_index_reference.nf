process p01_index_reference {
    tag "p01: bwa index on $reference"   
    publishDir "$params.outdir/p01_index_reference", mode: 'copy', overwrite: true

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
