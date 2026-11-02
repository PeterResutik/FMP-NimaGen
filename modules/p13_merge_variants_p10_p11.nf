process p13_merge_variants_p10_p11 {
    tag "p13: processing and merging variants on $sample_id"
    publishDir "$params.outdir/p13_merged_variants_xlsx", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(tssv_file), path(report_file), path(sc_file), path(sast_file), path(html_file)
    path reference
    path python_script_process_fdstools_sast
    path python_script_merge_fdstools_mutect2
    path python_modules // notation.py, frames.py, iupac.py: imported by the p13 scripts, inputs so -resume sees changes

    output:
    path("${sample_id}_${params.frame}_frame.xlsx"), emit: merged_variants_ch

    script:
    """
    if [ -s "${sast_file}" ]; then

        python $python_script_process_fdstools_sast \
            ${sast_file} \
            ${sample_id}_fdstools_processed.txt \
            --min_vf $params.min_vf_FDS --depth $params.depth --lh_thresh $params.lh_thresh --marker_map $params.fdstools_library \
            --tssv ${tssv_file} --reference $reference --frame $params.frame

        python $python_script_merge_fdstools_mutect2 \
            ${sample_id}_fdstools_processed.txt \
            ${sample_id}_${params.frame}_frame.xlsx \
            --min_vf $params.min_vf_FDS --marker_map $params.fdstools_library --variant_note $params.variant_note



    else
        echo "Skipping ${sample_id}: the FDSTOOLS output is empty." >&2
        python -c "import pandas as pd; pd.DataFrame([['This sample had no valid variant data']]).to_excel('${sample_id}_${params.frame}_frame.xlsx', index=False, header=False)"
    fi
    """
}
