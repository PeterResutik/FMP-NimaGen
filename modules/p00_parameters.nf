process p00_parameters{
    publishDir "$params.outdir",  mode: 'copy', overwrite: true

    input:
    val(logs)

    output:
    path "p00_parameters.txt"

    script:
    """
    echo '$logs' > p00_parameters.txt
    """
}
