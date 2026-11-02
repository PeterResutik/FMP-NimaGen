process p01b_prepare_humans_index {
    tag "p01b: prepare BWA index for humans reference (runs once)"

    input:
    path humans_index_dir
    val humans_base

    output:
    // the file name, once the index is ready: a value, so rtn's cache does not depend on this task's work folder
    val humans_base, emit: humans_index_ready

    script:
    """
    # The index is complete once .sa exists, the last file bwa index writes; an interrupted
    # unpacking or indexing leaves no complete-looking file behind and is redone
    if [[ ! -f "${humans_index_dir}/${humans_base}.sa" ]]; then
        echo "BWA index for human reference not found. Preparing it now..."
        if [[ ! -f "${humans_index_dir}/${humans_base}" ]]; then
            bunzip2 -c "${humans_index_dir}/${humans_base}.bz2" > "${humans_index_dir}/${humans_base}.tmp"
            mv "${humans_index_dir}/${humans_base}.tmp" "${humans_index_dir}/${humans_base}"
        fi
        bwa index "${humans_index_dir}/${humans_base}"
    else
        echo "BWA index found for human reference."
    fi
    # rtn reads the genomes through the .fai: built here once, not by the first parallel rtn runs
    if [[ ! -f "${humans_index_dir}/${humans_base}.fai" ]]; then
        samtools faidx "${humans_index_dir}/${humans_base}"
    fi
    """
}
