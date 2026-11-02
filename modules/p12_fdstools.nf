process p12_fdstools {
    tag "p12: fdstools on $sample_id"
    cpus { Math.min(params.num_threads as int, Runtime.runtime.availableProcessors()) }  // never more cores than the computer has
    publishDir { "$params.outdir/p12_fdstools/${sample_id}" }, mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(rtn_fastq)
    path fdstools_library

    output:
    tuple  val(sample_id), path("${sample_id}.tssv.csv"), path("${sample_id}.report.txt"), path("${sample_id}.sc.csv"), path("${sample_id}.sast.csv"), path("${sample_id}.html")

    script:
    // --baseQ is the lowest quality TSSV trusts, passed as the highest it distrusts (30 becomes 29);
    // 0 or 1 trusts every base: TSSV runs without --distrust-base-quality, which takes only 1 and up
    def distrust = (params.baseQ as int) > 1 ? "--distrust-base-quality ${(params.baseQ as int) - 1}" : ""
    """
    if [ -s "${rtn_fastq}" ]; then

    
        # bases below --baseQ are distrusted (--baseQ 0: none): a rarer sequence that differs from a more
        # abundant one of the same length only at such bases is counted for the more abundant one
        fdstools tssv $fdstools_library ${rtn_fastq} ${sample_id}.tssv.csv --minimum $params.minimum --num-threads ${task.cpus} --report ${sample_id}.report.txt \
            ${distrust}
        fdstools seqconvert allelename ${sample_id}.tssv.csv ${sample_id}.sc.csv --library $params.fdstools_library
        # samplestats names the alleles itself from the TSSV output: FDSTOOLS 2.2 cannot read back
        # N3107DEL, the name seqconvert writes for the rCRS placeholder at 3107
        fdstools samplestats ${sample_id}.tssv.csv ${sample_id}.sast.csv --library $params.fdstools_library -F allelename
        fdstools vis --min-abs $params.min_abs --min-pct-of-max $params.min_pct_of_max --min-pct-of-sum $params.min_pct_of_sum --allele-min-abs $params.allele_min_abs --allele-min-pct-of-max $params.allele_min_pct_of_max --allele-min-pct-of-sum $params.allele_min_pct_of_sum sample ${sample_id}.sast.csv ${sample_id}.html
    

    else
        echo "SKIPPED: ${rtn_fastq} is empty or missing." >&2
        touch ${sample_id}.tssv.csv ${sample_id}.report.txt ${sample_id}.sc.csv ${sample_id}.sast.csv ${sample_id}.html
    fi
    """
}
