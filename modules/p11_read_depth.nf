process p11_read_depth {
    tag "p11: read depth on $sample_id"
    publishDir { "${params.outdir}/p11_read_depth/${sample_id}" }, mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(read_depth_txt), path(read_depth_txt_numts)
    path python_script_generate_read_depth_plot

    output:
    path(read_depth_txt)
    path(read_depth_txt_numts)
    path("*read_depth_plot.png")

    script:
    """
    # Remove any pre-existing plot from a rerun (ensures reproducibility)
    rm -f ${sample_id}_read_depth_plot.png

    # Generate static read depth plot
    python $python_script_generate_read_depth_plot $read_depth_txt $read_depth_txt_numts ${sample_id}_read_depth_plot.png
    """
}
