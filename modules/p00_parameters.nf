process p00_parameters{
    publishDir "$params.outdir",  mode: 'copy', overwrite: true

    input:
    val(logs)

    output:
    path "p00_parameters.txt"

    // written by Nextflow itself, not by a shell command, so that the text may hold any character (a ' in a path)
    exec:
    task.workDir.resolve("p00_parameters.txt").text = logs + "\n"
}
