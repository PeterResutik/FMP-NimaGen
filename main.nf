nextflow.enable.dsl=2

params.reads = "$baseDir/raw_data/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/resources/rCRS/rCRS_NimaGen.fasta"

// flash
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 140 // default in FLASH is 65
params.max_mismatch_density = 0.25 //default in FLASH is 0.25

// params.multiqc = "$baseDir/multiqc"
params.publish_dir_mode = "symlink"
params.outdir = "results"
params.publish_bams = true // false: BAM/FASTQ intermediates are not copied to outdir (e.g. simulations); VCFs, FDSTOOLS files, QC and the merged Excel still are

params.adapter = 'ATCATAACAAAAAATTTCCACCAAA' // chrM:278-302, just before the poly-C tract at 303-315: R1 is cut after it (p02)

params.left_primers = "$baseDir/resources/primers/left_primers.fasta"
// params.left_primers_rc = "$baseDir/resources/primers/left_primers_rc.fasta"
params.right_primers_rc = "$baseDir/resources/primers/right_primers_rc.fasta"
params.amplicon_middle_positions = "$baseDir/resources/amplicon_bed/amplicons_bed.txt"
// params.force_sites = "$baseDir/force_sites.vcf.gz"

// cutadapt
params.quality_cutoff = 25
params.minimum_length = 60
params.maximum_length = 600 // longer than any merged 2x226 read, so products spanning two amplicons are kept
params.primer_min_overlap = 15 // bases of a primer that must be in the read when it is incomplete at the read end

// bwa mem in p07 (primer-trimmed merged reads, used by Mutect2)
params.clipping_penalty = "100,100" // bwa mem -L; bwa default is 5,5

params.humans_index_dir = "$baseDir/resources/rtn_files/humans"
params.humans_index_base = "humans_NimaGen_2026-10.fa" // built by resources/rtn_files/humans/build/ from public sources

params.numts_index_dir = "$baseDir/resources/rtn_files/numts"
params.numts_index_base = "Calabrese_Dayama_Smart_Numts_modified.fa"
params.rtn_positions = "$baseDir/resources/rtn_files/positions/recurrent_positions.tsv" // rtn sees these positions at base quality 2; none = off
params.rtn_memory = "4 GB" // memory per rtn task (p09), which holds the humans index (about 2 GB): Nextflow starts no more tasks than the computer's memory holds

// params.mtdna_database = "$baseDir/HelixMTdb_20200327_short.vcf.gz"
// params.fdstools_library = "$baseDir/resources/fdstools/mtNG_lib2_v211-flank.txt" 
params.fdstools_library = "$baseDir/resources/fdstools/mtNG_library_file.txt"

// rtn
params.mapQ = 30
params.skip_numt_filter = false // true skips rtn, e.g. for simulated reads, which contain no NUMTs; the MAPQ filter still applies

// fdstools
params.minimum = 2
params.num_threads = 6 // threads of FDSTOOLS tssv (p11), which reserves as many cores (at most the computer's); no faster beyond 6
params.min_abs = 2 
params.min_pct_of_max = 0 
params.min_pct_of_sum = 3 
params.allele_min_abs = 2 
params.allele_min_pct_of_max = 0 
params.allele_min_pct_of_sum = 3
 
// mutect2
params.baseQ = 30
params.callable_depth = 6
params.initial_tumor_lod = 0
params.tumor_lod_to_emit = 0
params.native_pair_hmm_threads = 1 // threads of Mutect2's pair HMM: more save almost no time on mtDNA
params.max_reads_per_alignment_start = 0
params.mutect2_memory = "4 GB" // memory per Mutect2 task (p12); Java may use 80% of it
params.force_alleles = "$baseDir/resources/mutect2/force_alleles.vcf" // alleles Mutect2 always evaluates (--alleles)
params.min_reads_per_strand = 0 // disables StrictStrandBiasFilter — this pipeline's RC-PCR + FLASH-merged reads are always single-orientation, so any nonzero value guarantees FILTER=strict_strand on every call regardless of true variant quality

params.python_script_remove_scb = "$baseDir/resources/scripts/remove_soft_clipped_bases_improved.py"
params.python_script_generate_read_depth_plot = "$baseDir/resources/scripts/generate_read_depth_plot_improved.py"
params.python_script_process_fdstools_sast = "$baseDir/resources/scripts/process_fdstools_output_improved_better.py"
params.python_script_process_mutect2_vcfgz = "$baseDir/resources/scripts/process_mutect2_output_improved.py"
params.python_script_merge_fdstools_mutect2 = "$baseDir/resources/scripts/merge_fdstools_mutect2_improved.py"
params.python_script_rtn_qualities = "$baseDir/resources/scripts/rtn_qualities.py"

// merging
params.depth = 20
params.min_vf_MT2 = 2
params.min_vf_FDS = 2
params.lh_thresh = 10 // symmetric: floor=10%, ceiling=1-10%=90% — below floor not reported, floor-ceiling reported as LHP (lowercase), above ceiling reported as major (uppercase)
params.disagreement_average = "plain" // plain | depth_weighted: how the merge averages both callers' frequencies when they disagree on major vs minor
params.frame = "separate" // separate | shared: how 16180-16193 and 300-315 are written (resources/scripts/frames.py)
params.mutect2_disabled_regions = "none" // none | all | names separated by commas (57-60,300-315,16180-16193): complex regions where only Mutect2's major calls count
params.variant_note = true // false: the report has no variant_note column

    // rm -r "$baseDir/work"
    // rm -r "$baseDir/results"
    // rm .nextflow.*

process p00_pipeline_parameters{
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

process p01_index_reference_fasta {
    tag "p01: bwa index on $reference"   
    publishDir "$params.outdir/p01_index", mode: 'copy', overwrite: true

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

process p02_map_raw_fastq_p01 {
    tag "p02: bwa mem on $sample_id"
    publishDir { "$params.outdir/p02_mapped_w_scb_bam/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*.bam*', enabled: params.publish_bams.toString().toBoolean()

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
    # -O: only the whole sequence counts; by default 3 matching bases at the end of any read are cut off
    # --action=retain: the anchor stays in R1, only the bases after it go
    cutadapt -a ${params.adapter} -O ${params.adapter.size()} --action=retain -o ${reads[0]}  tmp.fastq.gz 

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
    publishDir { "$params.outdir/p03_mapped_wo_scb_bam/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*sorted.bam*', enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(sam_r1), path(sam_r2)
    path python_script_remove_scb
    path reference
    path index_files
    path left_primers
    path right_primers

    output:
    tuple val(sample_id), path("${sample_id}_R1_wo_scb.bam"), path("${sample_id}_R2_wo_scb.bam")
    tuple val(sample_id), path("${sample_id}_R1_wo_scb_sorted.bam"), path("${sample_id}_R2_wo_scb_sorted.bam"), path("${sample_id}_R1_wo_scb_sorted.bam.bai"), path("${sample_id}_R2_wo_scb_sorted.bam.bai")
    tuple val(sample_id), path("${sample_id}_R1_R2_wo_scb_sorted.bam"), path("${sample_id}_R1_R2_wo_scb_sorted.bam.bai")
    tuple val(sample_id), path("${sample_id}_R1_wo_scb.fastq"), path("${sample_id}_R2_wo_scb.fastq"), emit: p03_fastq_files_wo_scb_ch


    script:
    """ 
    # soft-clipped bases go, and so do the few bases a read reaches past a primer boundary (read-through)
    python $python_script_remove_scb ${sam_r1} ${sample_id}_R1_wo_scb.sam ${reference} ${left_primers} ${right_primers}
    python $python_script_remove_scb ${sam_r2} ${sample_id}_R2_wo_scb.sam ${reference} ${left_primers} ${right_primers}

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
    # -t 1: with more threads FLASH writes the merged reads in a different order on every run
    flash ${fastq_r1_wo_scb} ${fastq_r2_wo_scb} -m $params.min_overlap -M $params.max_overlap -x $params.max_mismatch_density -O -t 1 -o ${sample_id}_wo_scb_merged
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
    # X instead of ^ and \$ (non-internal): a primer may be incomplete at the read end, as where the first bases of the
    # oligo carry an error that bwa clips and p03 removes, but never inside the read; at least primer_min_overlap of its bases must be there
    tr '^\$' XX < $left_primers > left_primers_noninternal.fasta
    tr '^\$' XX < $right_primers > right_primers_noninternal.fasta
    cutadapt -g file:left_primers_noninternal.fasta -O $params.primer_min_overlap -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_wo_scb_merged_trimmed_left.fastq $fastq_wo_scb_merged
    cutadapt -a file:right_primers_noninternal.fasta -O $params.primer_min_overlap -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_wo_scb_merged_trimmed_left_right.fastq ${sample_id}_wo_scb_merged_trimmed_left.fastq
    # quality trimming after the primers are removed, so that it cannot cut into a primer
    cutadapt -q $params.quality_cutoff -m $params.minimum_length -M $params.maximum_length -o ${sample_id}_wo_scb_merged_trimmed.fastq ${sample_id}_wo_scb_merged_trimmed_left_right.fastq
    """
}
        
process p06_map_merged_bam_p01_p04 {
    tag "p06: bwa mem on $sample_id"
    // publishDir "$params.outdir/p06_mapped_wo_scb_merged_bam/${sample_id}", mode: 'copy', pattern: '*.bam*'

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
        
process p07_map_merged_trimmed_bam_p01_p05 {
    tag "p07: bwa mem on $sample_id"
    // publishDir "$params.outdir/p07_mapped_wo_scb_merged_trimmed_bam/${sample_id}", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged_trimmed)
    path reference
    path index_files

    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${fastq_wo_scb_merged_trimmed.baseName}.bam"), path("${fastq_wo_scb_merged_trimmed.baseName}.bam.bai"), path("${fastq_wo_scb_merged_trimmed.baseName}_read_depth.txt")

    script:
    """
    # read ends are amplicon ends after primer trimming: clip only clearly wrong ends
    bwa mem -L ${params.clipping_penalty} ${reference} ${fastq_wo_scb_merged_trimmed} | samtools view -Sb - | samtools sort -o ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam

    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam  -o ${fastq_wo_scb_merged_trimmed.baseName}.bam
    samtools index ${fastq_wo_scb_merged_trimmed.baseName}.bam

    samtools depth -a -b $amplicon_middle_positions ${fastq_wo_scb_merged_trimmed.baseName}.bam > ${fastq_wo_scb_merged_trimmed.baseName}_read_depth.txt
    """
}

process p08_filter_numts_merged_fastq_p06 {
    tag "p08: reads rtn kept on $sample_id"
    publishDir { "$params.outdir/p08_filtered_numts_bam_for_fdstoold/${sample_id}" }, mode: 'copy', overwrite: true, enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(bam_wo_scb_merged), path(bam_wo_scb_merged_index), path(kept_names)

    output:
    tuple val(sample_id), path("${bam_wo_scb_merged.baseName}.rtn.bam"), path("${bam_wo_scb_merged.baseName}.rtn.bam.bai"), path("${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq")

    script:
    // the NUMT decision comes from p09, where rtn judges the same reads without their primers: the primers
    // are the kit's oligos, and where they carry a rare rCRS base (750A, 8860A) every read would differ from almost every genome
    """
    samtools view -b -h -q $params.mapQ -N ${kept_names} ${bam_wo_scb_merged} > ${bam_wo_scb_merged.baseName}.rtn.bam
    samtools index ${bam_wo_scb_merged.baseName}.rtn.bam
    samtools fastq ${bam_wo_scb_merged.baseName}.rtn.bam > ${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq
    """
}

process p09_filter_numts_trimmed_merged_bam_p07 {
    tag "p09: rtn on $sample_id"
    memory params.rtn_memory
    publishDir { "$params.outdir/p09_filtered_numts_bam_for_mutect2/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*.bam*', enabled: params.publish_bams.toString().toBoolean()

    input:
    tuple val(sample_id), path(bam_wo_scb_merged_trimmed), path(bam_wo_scb_merged_trimmed_index), path(read_depth_txt)
    path amplicon_middle_positions
    path humans_index
    val humans_base
    path numts_index
    val numts_base
    val rtn_lowq
    path python_script_rtn_qualities


    output:
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam"), path("${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam.bai"), path(read_depth_txt), path("${bam_wo_scb_merged_trimmed.baseName}_read_depth_wo_NUMTs.txt"), emit: bams
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}_kept_names.txt"), emit: kept_names

    script:
    // rtn -i: indels do not count toward a read's distance, as when rtn picks the closest genome
    def base = bam_wo_scb_merged_trimmed.baseName
    def rtn = "rtn -i -h \"${humans_index}/${humans_base}\" -n \"${numts_index}/${numts_base}\" -b"
    // rtn_lowq: rtn evaluates a copy of the reads with base quality 2 at these positions; its decision
    // (MAPQ, tags) is kept and the reads get their own qualities back, so Mutect2 sees the original ones
    def numt_filter = params.skip_numt_filter.toString().toBoolean() ?
        "cp ${bam_wo_scb_merged_trimmed} ${base}.rtn.bam && samtools index ${base}.rtn.bam" :
        rtn_lowq == 'none' ? "${rtn} ${bam_wo_scb_merged_trimmed}" :
        """set -o pipefail
    samtools view -h ${bam_wo_scb_merged_trimmed} | python ${python_script_rtn_qualities} lower ${rtn_lowq} > ${base}.lowq.sam
    samtools view -b -o ${base}.lowq.bam ${base}.lowq.sam && samtools index ${base}.lowq.bam
    ${rtn} ${base}.lowq.bam
    samtools view -h ${base}.lowq.rtn.bam | python ${python_script_rtn_qualities} restore ${bam_wo_scb_merged_trimmed} > ${base}.restored.sam
    samtools view -b -o ${base}.rtn.bam ${base}.restored.sam && samtools index ${base}.rtn.bam
    rm -f ${base}.lowq.* ${base}.restored.sam"""
    """
    ${numt_filter}

    samtools view -b -h -q $params.mapQ ${bam_wo_scb_merged_trimmed.baseName}.rtn.bam > ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam
    samtools index ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam
    samtools view ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam | cut -f1 | sort -u > ${bam_wo_scb_merged_trimmed.baseName}_kept_names.txt
    samtools depth -a -b $amplicon_middle_positions ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam > ${bam_wo_scb_merged_trimmed.baseName}_read_depth_wo_NUMTs.txt

    mv ${bam_wo_scb_merged_trimmed.baseName}.rtn.bam ${bam_wo_scb_merged_trimmed.baseName}.rtn.unfiltered.bam
    """
}

process p10_quality_control_p09 {
    tag "p10: read depth on $sample_id"
    publishDir { "${params.outdir}/p10_quality_control/${sample_id}" }, mode: 'copy', overwrite: true
    
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

process p11_variant_calling_fdstools_sast_p08 {
    tag "p11: fdstools on $sample_id"
    cpus { Math.min(params.num_threads as int, Runtime.runtime.availableProcessors()) }  // never more cores than the computer has
    publishDir { "$params.outdir/p11_fdstools/${sample_id}" }, mode: 'copy', overwrite: true
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(rtn_fastq)
    path fdstools_library

    output:
    tuple  val(sample_id), path("${sample_id}.tssv.csv"), path("${sample_id}.report.txt"), path("${sample_id}.sc.csv"), path("${sample_id}.sast.csv"), path("${sample_id}.html")

    script:
    """
    if [ -s "${rtn_fastq}" ]; then

    
        fdstools tssv $fdstools_library ${rtn_fastq} ${sample_id}.tssv.csv --minimum $params.minimum --num-threads ${task.cpus} --report ${sample_id}.report.txt
        fdstools seqconvert allelename ${sample_id}.tssv.csv ${sample_id}.sc.csv --library $params.fdstools_library
        fdstools samplestats ${sample_id}.sc.csv ${sample_id}.sast.csv
        fdstools vis --min-abs $params.min_abs --min-pct-of-max $params.min_pct_of_max --min-pct-of-sum $params.min_pct_of_sum --allele-min-abs $params.allele_min_abs --allele-min-pct-of-max $params.allele_min_pct_of_max --allele-min-pct-of-sum $params.allele_min_pct_of_sum sample ${sample_id}.sast.csv ${sample_id}.html
    

    else
        echo "SKIPPED: ${rtn_fastq} is empty or missing." >&2
        touch ${sample_id}.tssv.csv ${sample_id}.report.txt ${sample_id}.sc.csv ${sample_id}.sast.csv ${sample_id}.html
    fi
    """
}

process p12_variant_calling_mutect2_vcfgz_p01_p09 {
    tag "p12: mutect2 on $sample_id"
    memory params.mutect2_memory
    publishDir { "$params.outdir/p12_mutect2/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*.vcf.gz*'
    publishDir { "$params.outdir/p12_mutect2/${sample_id}" }, mode: 'copy', overwrite: true, pattern: '*bamout.bam*', enabled: params.publish_bams.toString().toBoolean()
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(read_depth_txt), path(read_depth_txt_numts) // the coverages do not need to be passed to mutect2 process
    path reference
    path fasta_index
    path mutect2_index
    path force_alleles

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

    # --alleles needs a bgzipped, indexed VCF; --force_alleles may be plain or bgzipped
    bcftools view -Oz -o forced_alleles.vcf.gz ${force_alleles}
    tabix -p vcf forced_alleles.vcf.gz

    # --max-mnp-distance 0: one record per substituted position, as FDSTOOLS reports them
    # pileup detection: candidates also from the aligned reads, not only from assembly;
    # proper-pair check off (merged reads are unpaired), same base quality as calling
    gatk --java-options "-Xmx${avail_mem}m" \
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
        --max-mnp-distance 0 \
        --pileup-detection true \
        --pileup-detection-enable-indel-pileup-calling true \
        --pileup-detection-proper-pair-read-badness false \
        --pileup-detection-snp-basequality-filter ${params.baseQ} \
        --alleles forced_alleles.vcf.gz \
        --bam-output ${bam_file.baseName}.bamout.bam \
        --tmp-dir tmp_${sample_id} \
        -I ${bam_file} \
        -O raw.vcf.gz 

    samtools sort -o ${bam_file.baseName}_sorted.bamout.bam ${bam_file.baseName}.bamout.bam
    samtools index ${bam_file.baseName}_sorted.bamout.bam 

    gatk  --java-options "-Xmx${avail_mem}m" \
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

    # every record goes on: --min_vf_MT2 and the length floor apply in p13, after the
    # records of one change (the alleles of a length variant) are added up
    mv ${bam_file.baseName}.norm.vcf.gz ${bam_file.baseName}.vcf.gz

    tabix -f ${bam_file.baseName}.vcf.gz 
    """
}
        // --base-quality-score-threshold ${params.baseQ} \

process p13_merge_variants_p10_p11 {
    tag "p13: processing and merging variants on $sample_id"
    publishDir "$params.outdir/p13_merged_variants_xlsx", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_file_idx), path(tssv_file), path(report_file), path(sc_file), path(sast_file), path(html_file), path(mutect2_depth_txt)
    path reference
    path python_script_process_mutect2_vcfgz
    path python_script_process_fdstools_sast
    path python_script_merge_fdstools_mutect2
    path python_modules // notation.py, frames.py, iupac.py: imported by the p13 scripts, inputs so -resume sees changes

    output:
    path("${sample_id}_${params.frame}_frame.xlsx"), emit: merged_variants_ch

    script:
    def vcf_name = "${vcf_file}".replaceAll('.vcf.gz', '')

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
                --min_vf $params.min_vf_MT2 --lh_thresh $params.lh_thresh --frame $params.frame --min_reads $params.minimum
        else
            printf "MUTECT2\tvf_MT2\trd_MT2\tMBQ\n" > ${vcf_file.baseName}.filtered.empop.txt
        fi

        python $python_script_process_fdstools_sast \
            ${sast_file} \
            ${sample_id}_fdstools_processed.txt \
            --min_vf $params.min_vf_FDS --depth $params.depth --lh_thresh $params.lh_thresh --marker_map $params.fdstools_library \
            --tssv ${tssv_file} --reference $reference --frame $params.frame

        python $python_script_merge_fdstools_mutect2 \
            ${sample_id}_fdstools_processed.txt \
            ${vcf_file.baseName}.filtered.empop.txt \
            ${sample_id}_${params.frame}_frame.xlsx \
            --lh_thresh $params.lh_thresh --min_vf $params.min_vf_FDS \
            --mutect2_depth ${mutect2_depth_txt} --marker_map $params.fdstools_library --depth $params.depth \
            --disagreement_average $params.disagreement_average \
            --mutect2_disabled_regions $params.mutect2_disabled_regions --min_reads $params.minimum \
            --variant_note $params.variant_note



    else
        echo "Skipping ${sample_id}: one or more input files are empty." >&2
        python -c "import pandas as pd; pd.DataFrame([['This sample had no valid variant data']]).to_excel('${sample_id}_${params.frame}_frame.xlsx', index=False, header=False)"
    fi
    """
}

// A path inside the pipeline folder, written relative to it, so the parameter log stays short
def relativePath(path) {
    return path.toString().replace("${projectDir}/", "")
}

workflow {
    log_text = """\
         m t D N A - N i m a G e n  P I P E L I N E    
         ==========================================
         mtDNA reference genome           : ${relativePath(params.reference)}
         --reads                          : ${relativePath(params.reads)}
         
         MERGING (with FLASH)
         --min_overlap                    : $params.min_overlap # The minimum required overlap length between two reads to provide a confident overlap (default: 10bp) 
         --max_overlap                    : $params.max_overlap # Maximum overlap length expected in approximately 90% of read pairs. 
         --max_mismatch_density           : $params.max_mismatch_density # Maximum allowed ratio between the number of mismatched base pairs and the overlap length 

         TRIMMING (with CUTADAPT) 
         --adapter                        : $params.adapter # Anchor before the poly-C tract (278-302); R1 is cut after it (p02)
         --left_primers                   : ${relativePath(params.left_primers)} # Left primers (soft clips p03, trimming p05)
         --right_primers_rc               : ${relativePath(params.right_primers_rc)} # Right primers, reverse complement (p03, p05)
         --quality_cutoff                 : $params.quality_cutoff # Trim low-quality bases from the end of each read, after primer removal.
         --primer_min_overlap             : $params.primer_min_overlap # Bases of a primer that must be in the read when it is incomplete at the read end.
         --minimum_length                 : $params.minimum_length # Reads shorter than this are discarded
         --maximum_length                 : $params.maximum_length # Reads longer than this are discarded
         --discard-untrimmed              : enabled (hard coded) # Discard reads that do not contain an adapter/primer.

         MAPPING (with BWA MEM, p07)
         --clipping_penalty               : $params.clipping_penalty # Clipping penalty (-L) for the primer-trimmed merged reads Mutect2 runs on (bwa default: 5,5)
         --amplicon_middle_positions      : ${relativePath(params.amplicon_middle_positions)} # Middle position of each amplicon, for its read depth

         NUMTs REMOVAL (with RTN)
         --humans_index_dir               : ${relativePath(params.humans_index_dir)} # Folder of the human mitogenomes rtn compares reads with
         --humans_index_base              : $params.humans_index_base # Human mitogenomes, unpacked and indexed on the first run
         --numts_index_dir                : ${relativePath(params.numts_index_dir)} # Folder of the NUMT sequences rtn compares reads with
         --numts_index_base               : $params.numts_index_base # NUMT sequences (from rtn)
         --rtn_positions                  : ${relativePath(params.rtn_positions)} # positions rtn sees at quality 2 (none = off)
         --mapQ                           : $params.mapQ # Used to filter out reads assigned as NUMTs by RTN 
         --skip_numt_filter               : $params.skip_numt_filter # true skips rtn (simulated reads contain no NUMTs); the MAPQ filter still applies
         --rtn_memory                     : $params.rtn_memory # memory per rtn task, which holds the humans index (about 2 GB)

         VARIANT CALLING (with FDSTOOLS)
         (tssv)
         --minimum                        : $params.minimum # Reads a sequence (FDSTOOLS) or a base (MUTECT2, rCRS included) needs
         --num_threads                    : $params.num_threads # Threads of FDSTOOLS tssv, at most the cores of the computer

         VARIANT CALLING (with MUTECT2)
         --baseQ                          : $params.baseQ # Minimum base quality required to consider a base for calling, also in Mutect2 pileup detection
         --callable_depth                 : $params.callable_depth # Minimum depth to be considered callable for Mutect stats. Does not affect genotyping
         --min_reads_per_strand           : $params.min_reads_per_strand # Minimum alt reads required on both forward and reverse strands
         --force_alleles                  : ${relativePath(params.force_alleles)} # VCF of alleles Mutect2 always evaluates
         --native_pair_hmm_threads        : $params.native_pair_hmm_threads # threads of the Mutect2 pair HMM; more save almost no time on mtDNA
         --mutect2_memory                 : $params.mutect2_memory # memory per Mutect2 task; Java may use 80 percent of it

         POSTPROCESSING
         --depth                          : $params.depth # Read depth below which an amplicon is reported as LOW; both callers keep their calls there
         --min_vf_MT2                     : $params.min_vf_MT2 # Minor variant frequency threshold MUTECT2
         --min_vf_FDS                     : $params.min_vf_FDS # Minor variant frequency threshold FDSTOOLS
         --lh_thresh                      : $params.lh_thresh # Length variants: under it not reported, up to 100 minus it lowercase, above uppercase
         --disagreement_average           : $params.disagreement_average # plain or depth_weighted average when the callers disagree on major vs minor
         --frame                          : $params.frame # separate (run lengths) or shared (laid from the left) for 16180-16193 and 300-315
         --mutect2_disabled_regions       : $params.mutect2_disabled_regions # none, all or e.g. 300-315,16180-16193: Mutect2 minor calls left out there
         --variant_note                   : $params.variant_note # false: the report has no variant_note column
         --marker_map                     : ${relativePath(params.fdstools_library)} # Path to marker map file

         OUTPUT DIRECTORY   
         --outdir                         : ${params.outdir}
         --publish_bams                   : $params.publish_bams # false: BAM/FASTQ intermediates are not copied to outdir
         """

    log.info(log_text)

    assert params.reference : "Missing reference genome path"
    assert params.reads : "Missing input reads"
    assert file(params.humans_index_dir).exists() : "Humans index directory does not exist"
    assert params.rtn_positions.toString() == 'none' || file(params.rtn_positions).exists() : "rtn positions file does not exist"

    channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    humans_index_ch = channel.value(file(params.humans_index_dir))
    numts_index_ch  = channel.value(file(params.numts_index_dir))

    // positions rtn sees at base quality 2 (resources/rtn_files/positions/README.md), comma-separated
    def rtn_lowq = params.rtn_positions.toString() == 'none' ? 'none' :
        file(params.rtn_positions).readLines().drop(1).collect { line -> line.split('\t')[1] }.join(',')
    humans_base_ch = channel.value(params.humans_index_base)
    numts_base_ch  = channel.value(params.numts_index_base)

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
    p03_filter_softclipped_fastq_p01_p02(p02_map_raw_fastq_p01.out.p02_raw_sam_ch, params.python_script_remove_scb, params.reference, p01_index_ch, params.left_primers, params.right_primers_rc)
    // p04_convert_bam_2_fastq_p03(p03_filter_softclipped_sam_p02.out.p03_bam_files_wo_scb_ch)

    // ──────────────── MERGING & TRIMMING ────────────────────
    p04_merge_fastq_p03(p03_filter_softclipped_fastq_p01_p02.out.p03_fastq_files_wo_scb_ch)
    p05_trim_merged_fastq_p04(p04_merge_fastq_p03.out, params.left_primers, params.right_primers_rc)

    // ────────────────── FINAL MAPPING ───────────────────────
    p06_map_merged_bam_p01_p04(p04_merge_fastq_p03.out, params.reference, p01_index_ch, params.amplicon_middle_positions)
    p07_map_merged_trimmed_bam_p01_p05(p05_trim_merged_fastq_p04.out, params.reference, p01_index_ch, params.amplicon_middle_positions)

    // ────────────────── NUMTs FILTERING ─────────────────────
    // rtn runs once, on the primer-trimmed reads (p09); p08 keeps the same reads for FDSTOOLS
    // rtn takes the humans folder itself and the file name from p01b once the index is ready
    p09_filter_numts_trimmed_merged_bam_p07(p07_map_merged_trimmed_bam_p01_p05.out, params.amplicon_middle_positions,  humans_index_ch, humans_index_ready_ch, numts_index_ch, numts_base_ch,
                                            rtn_lowq, params.python_script_rtn_qualities)
    p08_filter_numts_merged_fastq_p06(p06_map_merged_bam_p01_p04.out.join(p09_filter_numts_trimmed_merged_bam_p07.out.kept_names, by: 0))
    
    // ────────────────── QUALITY CONTROL ─────────────────────
    p10_quality_control_p09(p09_filter_numts_trimmed_merged_bam_p07.out.bams, params.python_script_generate_read_depth_plot)
    

    // ────────────────── VARIANT CALLING ─────────────────────
    p10_fdstools_ch = p11_variant_calling_fdstools_sast_p08(p08_filter_numts_merged_fastq_p06.out, params.fdstools_library)
    p12_variant_calling_mutect2_vcfgz_p01_p09(p09_filter_numts_trimmed_merged_bam_p07.out.bams, params.reference, p01_index_ch, p01_index_mutect2_ch, params.force_alleles)
    p11_mutect2_ch = p12_variant_calling_mutect2_vcfgz_p01_p09.out.mutect2_ch 


    // ────────────── PROCESS & MERGE VARIANTS ────────────────
    // Depth at each amplicon's middle position in the BAM Mutect2 runs on, for Mutect2's LOW
    p09_depth_ch = p09_filter_numts_trimmed_merged_bam_p07.out.bams.map { sid, _bam, _bai, _rd, rd_wo_numts -> tuple(sid, rd_wo_numts) }
    p10_p11_final_inputs = p11_mutect2_ch.join(p10_fdstools_ch, by: 0).join(p09_depth_ch, by: 0)
    p13_merge_variants_p10_p11(p10_p11_final_inputs, params.reference, params.python_script_process_mutect2_vcfgz, params.python_script_process_fdstools_sast, params.python_script_merge_fdstools_mutect2,
                               files("$baseDir/resources/scripts/{notation,frames,iupac}.py"))

}
