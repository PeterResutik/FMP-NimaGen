nextflow.enable.dsl=2

params.reads = "$baseDir/raw_data/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/resources/rCRS/rCRS_NimaGen.fasta"

// flash
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 140 // default in FLASH is 65
params.max_mismatch_density = 0.25 //default in FLASH is 0.25

// params.multiqc = "$baseDir/multiqc"
params.publish_dir_mode = "symlink"
params.outdir = "results_fixed"

params.adapter = 'ATCATAACAAAAAATTTCCACCAAA'

params.left_primers = "$baseDir/resources/primers/left_primers.fasta"
// params.left_primers_rc = "$baseDir/resources/primers/left_primers_rc.fasta"
params.right_primers_rc = "$baseDir/resources/primers/right_primers_rc.fasta"
params.amplicon_middle_positions = "$baseDir/resources/amplicon_bed/amplicons_bed.txt"
// params.force_sites = "$baseDir/force_sites.vcf.gz"

// cutadapt
params.quality_cutoff = 25
params.minimum_length = 60
params.maximum_length = 300

params.humans_index_dir = "$baseDir/resources/rtn_files/humans"
params.humans_index_base = "humans.fa"

params.numts_index_dir = "$baseDir/resources/rtn_files/numts"
params.numts_index_base = "Calabrese_Dayama_Smart_Numts_modified.fa"

// params.mtdna_database = "$baseDir/HelixMTdb_20200327_short.vcf.gz"
// params.fdstools_library = "$baseDir/resources/fdstools/mtNG_lib2_v211-flank.txt"
params.fdstools_library = "$baseDir/resources/fdstools/mtNG_library_file.txt"

// rtn
params.mapQ = 30

// fdstools
params.minimum = 2
params.num_threads = 6
params.min_reads_filt = 2
params.min_abs = 2
params.min_pct_of_max = 0
params.min_pct_of_sum = 3
params.allele_min_abs = 2
params.allele_min_pct_of_max = 0
params.allele_min_pct_of_sum = 3

// mutect2
params.detection_limit = 0.08
params.baseQ = 30
params.callable_depth = 10 // stats-only (GATK: "Does not affect genotyping") — aligned with params.depth, the same coverage-adequacy threshold FDSTOOLS' LOW flag uses, and with GATK's own default (previously left at 6 for no documented reason)
params.initial_tumor_lod = 0
params.tumor_lod_to_emit = 0
params.native_pair_hmm_threads = 4
params.max_reads_per_alignment_start = 0
params.min_reads_per_strand = 0 // disables StrictStrandBiasFilter — this pipeline's RC-PCR + FLASH-merged reads are always single-orientation, so any nonzero value guarantees FILTER=strict_strand on every call regardless of true variant quality
params.candidate_min_bq = 20 // base-quality floor for the bcftools pileup that nominates candidate sites for Mutect2 --alleles (p11). Only affects which sites get *proposed*; Mutect2 re-evaluates each with its own model at --min-base-quality-score, and unsupported ones are dropped by the AD[alt]>0 filter, so this can stay permissive.
params.min_alt_mbq = 20 // minimum median base quality of the ALT allele (p11). Force-calling bypasses Mutect2's own base-quality filtering, so without this a position where nearly every base is unusable can be reported as a confident heteroplasmy — chrM:297 in s25-04643-E2 (median BQ 13 on both alleles) came back as a 47% call. Genuine calls across all validation samples sit at 27-39, so 20 removes only that artifact.
params.max_mnp_distance = 0 // default is 1: Mutect2 merges phased substitutions within this distance into one combined MNP record instead of reporting each position separately. That bundling is an accident of assembly, not biology — validated on real samples (chrM:16180-16193, chrM:303-315) to fragment adjacent-position frequencies inconsistently (one position benefiting from two merged records, its neighbor stuck with only one) with zero regressions elsewhere in the genome across all tested samples. 0 = always report substitutions per-position, matching EMPOP's own per-position nomenclature.

params.python_script_remove_scb = "$baseDir/resources/scripts/remove_soft_clipped_bases_improved.py"
params.python_script_generate_read_depth_plot = "$baseDir/resources/scripts/generate_read_depth_plot_improved.py"
params.python_script_process_fdstools_sast = "$baseDir/resources/scripts/process_fdstools_output_improved_better.py"
params.python_script_process_mutect2_vcfgz = "$baseDir/resources/scripts/process_mutect2_output_improved.py"
params.python_script_merge_fdstools_mutect2 = "$baseDir/resources/scripts/merge_fdstools_mutect2_improved.py"
// Imported (not directly invoked) by process_fdstools_output_improved_better.py
// and process_mutect2_output_improved.py at runtime - must be staged alongside
// them in every process that runs either script, or the import fails.
params.python_script_repeat_regions = "$baseDir/resources/scripts/repeat_regions.py"
params.python_script_call_repeat_regions = "$baseDir/resources/scripts/call_repeat_regions.py"
params.python_script_read_evidence = "$baseDir/resources/scripts/read_evidence.py"

// merging
params.depth = 10
params.min_vf_MT2 = 5
params.min_vf_FDS = 5
params.lh_thresh = 10 // symmetric: floor=10%, ceiling=1-10%=90% — below floor not reported, floor-ceiling reported as LHP (lowercase), above ceiling reported as major (uppercase)
params.homopolymer_mutect2_reporting = "true_mutect2" // p13: how Mutect2's own calls inside the boundary-run region (leading run + different-base extension run, chrM:16180-16193 - the only one either mode touches, same restriction as the FDSTOOLS-side homopolymer_reporting above; every other reference-derived boundary-run region genome-wide is left as Mutect2 reported it) are handled. "bam_override" replaces Mutect2's rows there with direct BAM read counts (apply_homopolymer_region_overrides). "disabled" drops Mutect2's indel/length-axis rows there with no replacement, keeping ordinary substitution calls like T16189C untouched (disable_homopolymer_length_calls); merge_fdstools_mutect2_improved.py marks called_by_MUTECT2 as DISABLED for these. "true_mutect2" (default) leaves Mutect2's own calls as is. See process_mutect2_output_improved.py.
params.homopolymer_reporting = "boundary_run" // p13, FDSTOOLS side only, and restricted to just chrM 16180 through the end of its poly-C tract (16189 T>C's own consequence) - the only region any of this has been validated against; every other homopolymer run genome-wide (310, the other 8 boundary-run regions find_boundary_run_regions detects, any other plain 4+ run) is passed through exactly as FDSTOOLS reported it, for both values below. Chooses which underlying model reports homopolymer run lengths within that region: "shared_frame" bakes everything into one borrowed reference frame - not used by default, known to mis-report samples with real subpopulations at the run boundary (see report_boundary_run's docstring). "boundary_run" (default) splits the region's leading run and different-base extension run (chrM:16180-16183 / 16184-16193) into two independent axes. Either way, every row in the region also gets the single most common molecule's own genotype appended to its variant_note as a bracketed tag (e.g. "[dominant molecule (30.56%): A16182C, A16183C, T16189C, -16193.1C, -16193.2C]"), per W. Parson's reading of the ISFG/EMPOP convention (pers. comm. 2026-09-21) - the full distribution stays exactly as it always was, this is layered on top, not a replacement (see process_fdstools_output_improved_better.py).

    // rm -r "$baseDir/work"
    // rm -r "$baseDir/results"
    // rm .nextflow.*

log_text = """\
         m t D N A - N i m a G e n  P I P E L I N E
         ==========================================
         mtDNA reference genome           : ${params.reference}
         reads                            : ${params.reads}

         MERGING (with FLASH)
         --min_overlap                    : $params.min_overlap # The minimum required overlap length between two reads to provide a confident overlap (default: 10bp)
         --max_overlap                    : $params.max_overlap # Maximum overlap length expected in approximately 90% of read pairs.
         --max_mismatch_density           : $params.max_mismatch_density # Maximum allowed ratio between the number of mismatched base pairs and the overlap length

         TRIMMING (with CUTADAPT)
         --quality-cutoff                 : $params.quality_cutoff # Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal.
         --minimum-length                 : $params.minimum_length # Discard reads shorter than LEN. Default: 0
         --maximum-length                 : $params.maximum_length # Discard reads longer than LEN. Default: no limit
         --discard-untrimmed              : enabled (hard coded) # Discard reads that do not contain an adapter/primer.

         NUMTs REMOVAL (with RTN)
         --mapQ                           : $params.mapQ # Used to filter out reads assigned as NUMTs by RTN

         VARIANT CALLING (with FDSTOOLS)
         (tssv)
         --minimum                        : $params.minimum # report only sequences with this minimum number of reads (default: 2)
         --num_threads                    : $params.num_threads # number of worker threads to use (default: 1)
         (samplestats)
         --min_reads_filt                 : $params.min_reads_filt # the minimum number of reads (default: 1)

         VARIANT CALLING (with MUTECT2)
         --baseQ                          : $params.baseQ # Minimum base quality required to consider a base for calling
         --callable_depth                 : $params.callable_depth # Minimum depth to be considered callable for Mutect stats. Does not affect genotyping
         --min_reads_per_strand           : $params.min_reads_per_strand # Minimum alt reads required on both forward and reverse strands

         POSTPROCESSING
         --depth                          : $params.depth # Read depth threshold for low coverage
         --min_vf_MT2                     : $params.min_vf_MT2 # Minor variant frequency threshold MUTECT2
         --min_vf_FDS                     : $params.min_vf_FDS # Minor variant frequency threshold FDSTOOLS
         --lh_thresh                      : $params.lh_thresh # Symmetric length-heteroplasmy threshold: below it, not reported; between it and (1-it), reported as LHP (lowercase); above (1-it), reported as major (uppercase)
         --marker_map                     : $params.fdstools_library # Path to marker map file
         --homopolymer_reporting          : $params.homopolymer_reporting # How homopolymer run lengths are reported (FDSTOOLS side): shared_frame / boundary_run (default); dominant-molecule tag always included in variant_note
         --homopolymer_mutect2_reporting  : $params.homopolymer_mutect2_reporting # How Mutect2 calls in boundary-run regions are handled: bam_override / disabled / true_mutect2 (default)

         OUTPUT DIRECTORY
         outdir                           : ${params.outdir}
         """

log.info(log_text)

assert params.reference, "Missing reference genome path"
assert params.reads, "Missing input reads"
assert file(params.humans_index_dir).exists(), "Humans index directory does not exist"

process p00_pipeline_parameters{
    publishDir "$params.outdir",  mode: 'copy'

    input:
    val(logs)

    output:
    path "p00_parameters.txt"

    script:
    """
    echo '$logs' > p00_parameters.txt
    """
}

process p01_index_reference_fasta {
    tag "p01: bwa index on $reference"
    publishDir "$params.outdir/p01_index", mode: 'copy'

    input:
    path reference

    output:
    path("${reference}.*"), emit: reference
    path("${reference.baseName}.dict"), emit: reference_dict

    script:
    """
    bwa index $reference

    samtools faidx $reference
    samtools dict $reference -o ${reference.baseName}.dict
    """
}

process p01b_prepare_humans_index {
    tag "p01b: prepare BWA index for humans reference (runs once)"

    input:
    path humans_index_dir
    val humans_base

    output:
    path humans_index_dir, emit: humans_index_ready

    script:
    """
    # Check if BWA index for humans reference exists (one of the .amb/.bwt/.ann/.pac/.sa files)
    if [[ ! -f "${humans_index_dir}/${humans_base}.amb" ]]; then
        echo "BWA index for human reference not found. Preparing it now..."
        bunzip2 -k "${humans_index_dir}/humans.fa.bz2"
        bwa index "${humans_index_dir}/${humans_base}"
    else
        echo "BWA index found for human reference."
    fi
    """
}

process p02_map_raw_fastq_p01 {
    tag "p02: bwa mem on $sample_id"
    publishDir "$params.outdir/p02_mapped_w_scb_bam/${sample_id}", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(reads)
    path reference
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_R1.sam"), path("${sample_id}_R2.sam"), emit: p02_raw_sam_ch
    tuple path("${sample_id}_R1.bam"), path("${sample_id}_R2.bam"), path("${sample_id}_R1.bam.bai"), path("${sample_id}_R2.bam.bai")
    tuple path("${sample_id}_R1_R2.bam"), path("${sample_id}_R1_R2.bam.bai")

    script:
    """
    mv ${reads[0]} tmp.fastq.gz
    cutadapt -a ${params.adapter} -o ${reads[0]}  tmp.fastq.gz

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


process p03_filter_softclipped_fastq_p01_p02 {
    tag "p03: removing scb from $sample_id"
    publishDir "$params.outdir/p03_mapped_wo_scb_bam/${sample_id}", mode: 'copy', pattern: '*sorted.bam*'

    input:
    tuple val(sample_id), path(sam_r1), path(sam_r2)
    path python_script_remove_scb
    path reference
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_R1_wo_scb.bam"), path("${sample_id}_R2_wo_scb.bam")
    tuple val(sample_id), path("${sample_id}_R1_wo_scb_sorted.bam"), path("${sample_id}_R2_wo_scb_sorted.bam"), path("${sample_id}_R1_wo_scb_sorted.bam.bai"), path("${sample_id}_R2_wo_scb_sorted.bam.bai")
    tuple val(sample_id), path("${sample_id}_R1_R2_wo_scb_sorted.bam"), path("${sample_id}_R1_R2_wo_scb_sorted.bam.bai")
    tuple val(sample_id), path("${sample_id}_R1_wo_scb.fastq"), path("${sample_id}_R2_wo_scb.fastq"), emit: p03_fastq_files_wo_scb_ch


    script:
    """
    python $python_script_remove_scb ${sam_r1} ${sample_id}_R1_wo_scb.sam
    python $python_script_remove_scb ${sam_r2} ${sample_id}_R2_wo_scb.sam

    samtools view -Sb ${sample_id}_R1_wo_scb.sam > ${sample_id}_R1_wo_scb.bam
    samtools view -Sb ${sample_id}_R2_wo_scb.sam > ${sample_id}_R2_wo_scb.bam

    samtools fastq ${sample_id}_R1_wo_scb.bam > ${sample_id}_R1_wo_scb.fastq
    samtools fastq ${sample_id}_R2_wo_scb.bam > ${sample_id}_R2_wo_scb.fastq

    bwa mem ${reference} ${sample_id}_R1_wo_scb.fastq ${sample_id}_R2_wo_scb.fastq  > ${sample_id}_R1_R2_wo_scb.sam
    samtools view -bS ${sample_id}_R1_R2_wo_scb.sam | samtools sort -o ${sample_id}_R1_R2_wo_scb_sorted.bam
    samtools index ${sample_id}_R1_R2_wo_scb_sorted.bam

    samtools sort -o ${sample_id}_R1_wo_scb_sorted.bam ${sample_id}_R1_wo_scb.bam
    samtools sort -o ${sample_id}_R2_wo_scb_sorted.bam ${sample_id}_R2_wo_scb.bam
    samtools index ${sample_id}_R1_wo_scb_sorted.bam
    samtools index ${sample_id}_R2_wo_scb_sorted.bam
    """
}

process p04_merge_fastq_p03 {
    tag "p04: flash on $sample_id"
    // publishDir "$params.outdir/p04_wo_scb_merged_fastq/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_r1_wo_scb), path(fastq_r2_wo_scb)

    output:
    tuple val(sample_id), path("${sample_id}_wo_scb_merged.extendedFrags.fastq")

    script:
    """
    flash ${fastq_r1_wo_scb} ${fastq_r2_wo_scb} -m $params.min_overlap -M $params.max_overlap -x $params.max_mismatch_density -O -o ${sample_id}_wo_scb_merged
    """
}

process p05_trim_merged_fastq_p04 {
    tag "p05: cutadapt on $sample_id"
    // publishDir "$params.outdir/p05_wo_scb_merged_trimmed_fastq/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged)
    path left_primers
    path right_primers

    output:
    tuple val(sample_id), path("${sample_id}_wo_scb_merged_trimmed.fastq")

    script:
    """
    cutadapt -g file:$left_primers -q $params.quality_cutoff -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_wo_scb_merged_trimmed_left.fastq $fastq_wo_scb_merged
    cutadapt -a file:$right_primers -q $params.quality_cutoff -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_wo_scb_merged_trimmed_left_right.fastq ${sample_id}_wo_scb_merged_trimmed_left.fastq
    cp ${sample_id}_wo_scb_merged_trimmed_left_right.fastq ${sample_id}_wo_scb_merged_trimmed.fastq
    """
}

process p06_map_merged_trimmed_bam_p01_p05 {
    tag "p06: bwa mem on $sample_id"
    // publishDir "$params.outdir/p06_mapped_wo_scb_merged_trimmed_bam/${sample_id}", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged_trimmed)
    path reference
    path index_files

    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${fastq_wo_scb_merged_trimmed.baseName}.bam"), path("${fastq_wo_scb_merged_trimmed.baseName}.bam.bai"), path("${fastq_wo_scb_merged_trimmed.baseName}_read_depth.txt")

    script:
    """
    bwa mem ${reference} ${fastq_wo_scb_merged_trimmed} | samtools view -Sb - | samtools sort -o ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam

    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam  -o ${fastq_wo_scb_merged_trimmed.baseName}.bam
    samtools index ${fastq_wo_scb_merged_trimmed.baseName}.bam

    samtools depth -a -b $amplicon_middle_positions ${fastq_wo_scb_merged_trimmed.baseName}.bam > ${fastq_wo_scb_merged_trimmed.baseName}_read_depth.txt
    """
}

process p07_map_merged_bam_p01_p04 {
    tag "p07: bwa mem on $sample_id"
    // publishDir "$params.outdir/p07_mapped_wo_scb_merged_bam/${sample_id}", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged)
    path reference
    path index_files
    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${sample_id}_wo_scb_merged.bam"), path("${sample_id}_wo_scb_merged.bam.bai")

    script:
    """
    bwa mem ${reference} ${fastq_wo_scb_merged} | samtools view -Sb - | samtools sort -o ${sample_id}_wo_scb_sorted.bam

    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${sample_id}_wo_scb_sorted.bam  -o ${sample_id}_wo_scb_merged.bam
    samtools index ${sample_id}_wo_scb_merged.bam
    """
}

process p08_filter_numts_trimmed_merged_bam_p06 {
    tag "p08: rtn on $sample_id"
    publishDir "$params.outdir/p08_filtered_numts_bam_for_mutect2/${sample_id}", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(bam_wo_scb_merged_trimmed), path(bam_wo_scb_merged_trimmed_index), path(read_depth_txt)
    path amplicon_middle_positions
    path humans_index
    val humans_base
    path numts_index
    val numts_base


    output:
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam"), path("${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam.bai"), path(read_depth_txt), path("${bam_wo_scb_merged_trimmed.baseName}_read_depth_wo_NUMTs.txt"), emit: main_ch
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam"), path("${bam_wo_scb_merged_trimmed.baseName}.rtn.unfiltered.bam"), emit: numt_classification_ch

    script:
    """
    rtn -h "${humans_index}/${humans_base}" -n "${numts_index}/${numts_base}" -b $bam_wo_scb_merged_trimmed

    samtools view -b -h -q $params.mapQ ${bam_wo_scb_merged_trimmed.baseName}.rtn.bam > ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam
    samtools index ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam
    samtools depth -a -b $amplicon_middle_positions ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam > ${bam_wo_scb_merged_trimmed.baseName}_read_depth_wo_NUMTs.txt

    mv ${bam_wo_scb_merged_trimmed.baseName}.rtn.bam ${bam_wo_scb_merged_trimmed.baseName}.rtn.unfiltered.bam
    """
}

process p09_filter_numts_merged_fastq_p07 {
    tag "p09: rtn on $sample_id"
    publishDir "$params.outdir/p09_filtered_numts_bam_for_fdstoold/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_wo_scb_merged), path(bam_wo_scb_merged_index), path(trimmed_passed_rtn_bam), path(trimmed_unfiltered_rtn_bam)
    path amplicon_middle_positions
    path humans_index
    val humans_base
    path numts_index
    val numts_base

    output:
    tuple val(sample_id), path("${bam_wo_scb_merged.baseName}.rtn.bam"), path("${bam_wo_scb_merged.baseName}.rtn.bam.bai"), path("${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq")

    script:
    """
    rtn -h "${humans_index}/${humans_base}" -n "${numts_index}/${numts_base}" -b $bam_wo_scb_merged

    # NUMT/MAPQ classification decided here (on UNTRIMMED reads) is exposed
    # to primer-region sequencing errors cutadapt hasn't removed yet, which
    # degrades rtn's realignment confidence and zeroes MAPQ far more often
    # than it should - confirmed on HG00376: insertion-carrying reads in
    # the 16193 homopolymer are ~5.5x more likely to get MAPQ=0'd than the
    # region's reads overall (40.9% vs 7.5%), starving FDSTOOLS of real
    # length-heteroplasmy support that Mutect2's trimmed branch (p08) sees
    # fine. So: prefer p08's classification - computed on the SAME
    # underlying merged molecules, just after cutadapt has already removed
    # the error-prone primer bases - wherever a read has one. A read only
    # falls back to this run's own MAPQ decision when cutadapt discarded
    # it during trimming (p05's --discard-untrimmed) and p08 never saw it
    # at all, so there's no better-informed decision available for it.
    samtools view -q $params.mapQ ${bam_wo_scb_merged.baseName}.rtn.bam | cut -f1 | sort -u > own_passed_ids.txt
    samtools view $trimmed_unfiltered_rtn_bam | cut -f1 | sort -u > all_trimmed_ids.txt
    samtools view $trimmed_passed_rtn_bam | cut -f1 | sort -u > passed_trimmed_ids.txt
    comm -23 own_passed_ids.txt all_trimmed_ids.txt > fallback_ids.txt
    cat passed_trimmed_ids.txt fallback_ids.txt | sort -u > final_keep_ids.txt

    samtools view -h -b -N final_keep_ids.txt ${bam_wo_scb_merged.baseName}.rtn.bam > ${bam_wo_scb_merged.baseName}.rtn_tmp.bam
    samtools fastq ${bam_wo_scb_merged.baseName}.rtn_tmp.bam > ${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq
    """
}

process p10_quality_control_p08 {
    tag "p10: fastqc + read depth on $sample_id"
    publishDir "${params.outdir}/p10_quality_control/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(read_depth_txt), path(read_depth_txt_numts)
    path python_script_generate_read_depth_plot

    output:
    path "*.zip", emit: fastqc_ch
    path(read_depth_txt)
    path(read_depth_txt_numts)
    path("*read_depth_plot.png")

    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    """
    fastqc --threads ${task.cpus} --memory ${avail_mem} $bam_file -o .

    # Remove any pre-existing plot from a rerun (ensures reproducibility)
    rm -f ${sample_id}_read_depth_plot.png

    # Generate static read depth plot
    python $python_script_generate_read_depth_plot $read_depth_txt $read_depth_txt_numts ${sample_id}_read_depth_plot.png
    """
}

process p11_variant_calling_mutect2_vcfgz_p01_p08 {
    tag "p11: mutect2 on $sample_id"
    publishDir "$params.outdir/p11_mutect2/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(read_depth_txt), path(read_depth_txt_numts) // the coverages do not need to be passed to mutect2 process
    path reference
    path fasta_index
    path mutect2_index

    output:
    tuple  val(sample_id), path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), emit: mutect2_ch
    path("${bam_file.baseName}_sorted.bamout.bam")
    path("${bam_file.baseName}_sorted.bamout.bam.bai")


    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    mkdir tmp_${sample_id}

    ## Nominate candidate sites from a pileup-based caller before running
    ## Mutect2. Mutect2's local re-assembly can silently fail to *propose*
    ## an allele at all in this panel's densely tiled, high-depth amplicon
    ## regions - not a scoring/threshold issue but a de Bruijn graph one, so
    ## no LOD/sensitivity flag recovers it (verified: mitochondria-mode,
    ## genotype-germline-sites, force-active, pileup-detection, adaptive
    ## pruning off, kmer sweeps, shifted reference all still miss them).
    ## Empirically a 5bp shift in the assembly window flips whether a
    ## haplotype with 100% read support survives, so no parameter value
    ## fixes it generally. bcftools has no assembly step and sees these
    ## sites trivially; passing them via --alleles makes Mutect2 genotype
    ## them with its own model (it still runs normal discovery as well).
    ## Rescues real variants default Mutect2 drops entirely - e.g. a 99%
    ## 573.1C insertion and a 47% A297C heteroplasmy, both confirmed
    ## against raw reads and FDSTOOLS. See broadinstitute/gatk#9137.
    bcftools mpileup \
        -f ${reference} \
        -d 100000 -L 100000 -a AD \
        -Q ${params.candidate_min_bq} -q 0 -Ou \
        ${bam_file} \
        | bcftools call -mv -Ov \
        | bcftools norm -f ${reference} -m-any -Oz -o candidates.vcf.gz

    tabix -f candidates.vcf.gz

    gatk --java-options "-Xmx16G" \
        Mutect2 \
        -R ${reference} \
        -L 'chrM' \
        --min-base-quality-score ${params.baseQ} \
        --callable-depth $params.callable_depth \
        --linked-de-bruijn-graph true \
        --recover-all-dangling-branches true \
        --initial-tumor-lod "${params.initial_tumor_lod}" \
        --tumor-lod-to-emit "${params.tumor_lod_to_emit}"  \
        --native-pair-hmm-threads "${params.native_pair_hmm_threads}"  \
        --max-reads-per-alignment-start "${params.max_reads_per_alignment_start}" \
        --max-mnp-distance "${params.max_mnp_distance}" \
        --alleles candidates.vcf.gz \
        --bam-output ${bam_file.baseName}.bamout.bam \
        --tmp-dir tmp_${sample_id} \
        -I ${bam_file} \
        -O raw.vcf.gz

    samtools sort -o ${bam_file.baseName}_sorted.bamout.bam ${bam_file.baseName}.bamout.bam
    samtools index ${bam_file.baseName}_sorted.bamout.bam

    gatk  --java-options "-Xmx16G" \
        FilterMutectCalls \
        -R ${reference} \
        --min-reads-per-strand "${params.min_reads_per_strand}" \
        -V raw.vcf.gz \
        --tmp-dir . \
        -O ${bam_file.baseName}.vcf.gz

    bcftools norm \
        -m-any \
        -f ${reference} \
        -o ${bam_file.baseName}.norm.vcf.gz -Oz \
        ${bam_file.baseName}.vcf.gz

    detection_limit=\$(echo "${params.min_vf_MT2} / 100" | bc -l)

    ## AD[0:1]>=min_reads_filt drops force-called sites Mutect2 evaluated
    ## and found too little support for. The >0 case is required, not
    ## optional: for zero-support sites Mutect2 reports a degenerate
    ## prior-ish AF (0.5, 0.406 observed) despite AD=x,0, so an AF
    ## threshold alone lets them through. Raised from >0 to
    ## >=min_reads_filt (2) to also drop single-read "heteroplasmy" calls
    ## at low depth - the same protection FDSTOOLS already has via its own
    ## --min-reads-filt singleton removal (samplestats step), just missing
    ## on the Mutect2 side. Validated against all 15 samples: removes only
    ## HG01799's T4950C (AD=15,1, already flagged as a Mutect2-only false
    ## positive) and s24-12883a's two 1-total-read calls at 16189/16192 -
    ## nothing else in the validation set has alt support below 2 reads.
    ##
    ## MBQ[1] is the median base quality of the ALT allele, and guards the
    ## other failure mode of force-calling: it bypasses the base-quality
    ## filter Mutect2's own discovery applies, so a site where nearly every
    ## base is unusable can come back with a confident-looking AD/AF. In
    ## s25-04643-E2, chrM:297 has median BQ 13 on BOTH alleles (only 3 of 26
    ## alt reads reach Q30) - normal discovery correctly skipped it, but
    ## force-calling reported it as a 47% heteroplasmy with AD=23,22. An AD
    ## filter cannot catch this because that AD counts the low-quality reads.
    ## Raising the nomination -Q does not help either: the reference reads
    ## are equally bad, so filtering them makes the site look stronger.
    ## The threshold sits in a natural gap in the data - across all
    ## validation samples the only alt MBQ below 20 was that artifact (13),
    ## while genuine calls run 27-39, including real ones at 27-29 (e.g. a
    ## 99.8% homoplasmic G8557A on 454 reads) that a Q30 cut would discard.
    bcftools view \
    -i "FORMAT/AF>=\$detection_limit && FORMAT/AD[0:1]>=${params.min_reads_filt} && MBQ[1]>=${params.min_alt_mbq}" \\
    -o ${bam_file.baseName}.vcf.gz -Oz \
    ${bam_file.baseName}.norm.vcf.gz

    tabix -f ${bam_file.baseName}.vcf.gz
    """
}
        // --base-quality-score-threshold ${params.baseQ} \

process p12_variant_calling_fdstools_sast_p09 {
    tag "p12: fdstools on $sample_id"
    publishDir "$params.outdir/p12_fdstools/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(rtn_fastq)
    path fdstools_library

    output:
    tuple  val(sample_id), path("${sample_id}.tssv.csv"), path("${sample_id}.report.txt"), path("${sample_id}.sc.csv"), path("${sample_id}.sast.csv"), path("${sample_id}.html")

    script:
    """
    if [ -s "${rtn_fastq}" ]; then


        fdstools tssv $fdstools_library ${rtn_fastq} ${sample_id}.tssv.csv --minimum $params.minimum --num-threads $params.num_threads --report ${sample_id}.report.txt
        fdstools seqconvert allelename ${sample_id}.tssv.csv ${sample_id}.sc.csv --library $params.fdstools_library
        fdstools samplestats --min-reads-filt $params.min_reads_filt ${sample_id}.sc.csv ${sample_id}.sast.csv
        fdstools vis --min-abs $params.min_abs --min-pct-of-max $params.min_pct_of_max --min-pct-of-sum $params.min_pct_of_sum --allele-min-abs $params.allele_min_abs --allele-min-pct-of-max $params.allele_min_pct_of_max --allele-min-pct-of-sum $params.allele_min_pct_of_sum sample ${sample_id}.sast.csv ${sample_id}.html


    else
        echo "SKIPPED: ${rtn_fastq} is empty or missing." >&2
        touch ${sample_id}.tssv.csv ${sample_id}.report.txt ${sample_id}.sc.csv ${sample_id}.sast.csv ${sample_id}.html
    fi
    """
}

process p13_merge_variants_p11_p12 {
    tag "p13: processing and merging variants on $sample_id"
    publishDir "$params.outdir/p13_merged_variants_xlsx", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_file_idx), path(tssv_file), path(report_file), path(sc_file), path(sast_file), path(html_file), path(mutect2_bam), path(mutect2_bam_idx)
    path reference
    path python_script_process_mutect2_vcfgz
    path python_script_process_fdstools_sast
    path python_script_merge_fdstools_mutect2
    path python_script_repeat_regions
    path python_script_call_repeat_regions
    path python_script_read_evidence

    output:
    path("${sample_id}_merged_variants.xlsx"), emit: merged_variants_ch

    script:
    def vcf_name = "${vcf_file}".replaceAll('.vcf.gz', '')
    def bam_arg = "--bam ${mutect2_bam} --homopolymer_mutect2_reporting ${params.homopolymer_mutect2_reporting}"

    """
    if [ -s "${vcf_file}" ] && [ -s "${sast_file}" ]; then

        echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tMeanBaseQuality\tCoverage\tGT" \
            > ${vcf_file.baseName}.txt

        bcftools query -u \
            -f "${vcf_name}.bam\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%MBQ\t%AD\t%GT]\n" \
            ${vcf_file} | \
            awk -F'\\t' '(\$2 !~ /bla/)' \
            >> ${vcf_file.baseName}.txt

        ## Annotate SNVs and INDELs
        awk 'BEGIN {OFS="\\t"} {
            if (NR == 1) { print \$0, "Type"; next }
            if ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5)) { \$10="INDEL" }
            else { \$10="SNP" }
            print
        }' ${vcf_file.baseName}.txt > ${vcf_file.baseName}.filtered.txt

        mutect_rows=\$(awk 'NR>1{c++} END{print c+0}' ${vcf_file.baseName}.filtered.txt)

        if [ "\$mutect_rows" -gt 0 ]; then
            python $python_script_process_mutect2_vcfgz \
                ${vcf_file.baseName}.filtered.txt \
                ${vcf_file.baseName}.filtered.empop.txt \
                $reference \
                --min_vf $params.min_vf_MT2 --lh_thresh $params.lh_thresh \
                $bam_arg
        else
            printf "MUTECT2\tvf_MT2\trd_MT2\tMBQ\n" > ${vcf_file.baseName}.filtered.empop.txt
        fi

        python $python_script_process_fdstools_sast \
            ${sast_file} \
            ${sample_id}_fdstools_processed.txt \
            --reference $reference \
            --min_vf $params.min_vf_FDS --depth $params.depth --lh_thresh $params.lh_thresh --marker_map $params.fdstools_library \
            --homopolymer_reporting $params.homopolymer_reporting

        python $python_script_merge_fdstools_mutect2 \
            ${sample_id}_fdstools_processed.txt \
            ${vcf_file.baseName}.filtered.empop.txt \
            ${sample_id}_merged_variants.xlsx \
            --lh_thresh $params.lh_thresh --min_vf $params.min_vf_FDS \
            --homopolymer_mutect2_reporting $params.homopolymer_mutect2_reporting --reference $reference



    else
        echo "Skipping ${sample_id}: one or more input files are empty." >&2
        python -c "import pandas as pd; pd.DataFrame([['This sample had no valid variant data']]).to_excel('${sample_id}_merged_variants.xlsx', index=False, header=False)"
    fi
    """
}

workflow {
    // p13's own output gets deleted and freshly republished on every run,
    // regardless of -resume - reusing the same --outdir across two runs
    // that differ only in homopolymer_reporting/homopolymer_mutect2_
    // reporting (or any other param p13 depends on) can otherwise leave a
    // STALE file in place: publishDir doesn't reliably re-copy a cached
    // task's output when something already sits at that exact
    // destination path from an earlier run with different params, so the
    // merged xlsx can keep showing an old configuration's results even
    // though the console reports the (differently-configured) task as
    // correctly "Cached" (2026-09-23, user, after --outdir results_
    // boundary_run_new kept showing "disabled" output after switching
    // back to true_mutect2 with the same --outdir: "couldn't we do
    // something like delete the p13 folder in the pipeline for every
    // rerun?"). Only p13's own directory is cleared - every earlier
    // stage's output is untouched, and -resume still skips genuinely
    // unchanged upstream work exactly as before; this only forces a
    // fresh *publish* of whatever p13 already has cached or computes.
    def p13_outdir = file("${params.outdir}/p13_merged_variants_xlsx")
    if (p13_outdir.exists()) {
        p13_outdir.deleteDir()
    }

    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    humans_index_ch = Channel.value(file(params.humans_index_dir))
    numts_index_ch  = Channel.value(file(params.numts_index_dir))

    humans_base_ch = Channel.value(params.humans_index_base)
    numts_base_ch  = Channel.value(params.numts_index_base)

    // ────────────────── LOG PARAMETERS ──────────────────────
    p00_pipeline_parameters(log_text)

    // ──────────────── REFERENCE INDEXING ────────────────────
    p01_index_reference_fasta(params.reference)
    p01_index_ch = p01_index_reference_fasta.out.reference
    p01_index_mutect2_ch = p01_index_reference_fasta.out.reference_dict

    // ──────── HUMANS REFERENCE INDEXING (runs once) ─────────
    p01b_prepare_humans_index(humans_index_ch, humans_base_ch)
    humans_index_ready_ch = p01b_prepare_humans_index.out.humans_index_ready

    // ─────────── RAW READ MAPPING & FILTERING ───────────────
    p02_map_raw_fastq_p01(read_pairs_ch, params.reference, p01_index_ch)
    p03_filter_softclipped_fastq_p01_p02(p02_map_raw_fastq_p01.out.p02_raw_sam_ch, params.python_script_remove_scb, params.reference, p01_index_ch)
    // p04_convert_bam_2_fastq_p03(p03_filter_softclipped_sam_p02.out.p03_bam_files_wo_scb_ch)

    // ──────────────── MERGING & TRIMMING ────────────────────
    p04_merge_fastq_p03(p03_filter_softclipped_fastq_p01_p02.out.p03_fastq_files_wo_scb_ch)
    p05_trim_merged_fastq_p04(p04_merge_fastq_p03.out, params.left_primers, params.right_primers_rc)

    // ────────────────── FINAL MAPPING ───────────────────────
    p06_map_merged_trimmed_bam_p01_p05(p05_trim_merged_fastq_p04.out, params.reference, p01_index_ch, params.amplicon_middle_positions)
    p07_map_merged_bam_p01_p04(p04_merge_fastq_p03.out, params.reference, p01_index_ch, params.amplicon_middle_positions)

    // ────────────────── NUMTs FILTERING ─────────────────────
    // p08 (trimmed branch) runs first: p09 (untrimmed/FDSTOOLS branch)
    // borrows its NUMT/MAPQ classification instead of trusting its own
    // (see p09's script for why) - so p09 now depends on p08's output.
    p08_filter_numts_trimmed_merged_bam_p06(p06_map_merged_trimmed_bam_p01_p05.out, params.amplicon_middle_positions,  humans_index_ready_ch, humans_base_ch, numts_index_ch, numts_base_ch)
    p09_inputs = p07_map_merged_bam_p01_p04.out.join(p08_filter_numts_trimmed_merged_bam_p06.out.numt_classification_ch, by: 0)
    p09_filter_numts_merged_fastq_p07(p09_inputs, params.amplicon_middle_positions, humans_index_ready_ch, humans_base_ch, numts_index_ch, numts_base_ch)

    // ────────────────── QUALITY CONTROL ─────────────────────
    p10_quality_control_p08(p08_filter_numts_trimmed_merged_bam_p06.out.main_ch, params.python_script_generate_read_depth_plot)


    // ────────────────── VARIANT CALLING ─────────────────────
    p12_fdstools_ch = p12_variant_calling_fdstools_sast_p09(p09_filter_numts_merged_fastq_p07.out, params.fdstools_library)
    p11_variant_calling_mutect2_vcfgz_p01_p08(p08_filter_numts_trimmed_merged_bam_p06.out.main_ch, params.reference, p01_index_ch, p01_index_mutect2_ch)
    p11_mutect2_ch = p11_variant_calling_mutect2_vcfgz_p01_p08.out.mutect2_ch


    // ────────────── PROCESS & MERGE VARIANTS ────────────────
    // Same BAM Mutect2 itself was run on (p08's filtered/mapQ output) -
    // p13 needs it directly too, to read repeat-region loci straight
    // from the reads rather than trusting Mutect2's own VCF bookkeeping
    // there (see process_mutect2_output_improved.py).
    p08_bam_for_p13 = p08_filter_numts_trimmed_merged_bam_p06.out.main_ch.map { sid, bam, bai, rd1, rd2 -> tuple(sid, bam, bai) }
    p11_p12_final_inputs = p11_mutect2_ch.join(p12_fdstools_ch, by: 0).join(p08_bam_for_p13, by: 0)
    p13_merge_variants_p11_p12(p11_p12_final_inputs, params.reference, params.python_script_process_mutect2_vcfgz, params.python_script_process_fdstools_sast, params.python_script_merge_fdstools_mutect2, params.python_script_repeat_regions, params.python_script_call_repeat_regions, params.python_script_read_evidence)

}
