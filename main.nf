nextflow.enable.dsl=2

params.reads = "$baseDir/raw_data/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/resources/rCRS/rCRS_NimaGen.fasta"

// flash
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 140 // default in FLASH is 65
params.max_mismatch_density = 0.25 //default in FLASH is 0.25

// params.multiqc = "$baseDir/multiqc"
params.publish_dir_mode = "symlink"
params.outdir = "results_new"

params.adapter = 'ATCATAACAAAAAATTTCCACCAAA'

params.left_primers = "$baseDir/resources/primers/left_primers.fasta"
params.left_primers_rc = "$baseDir/resources/primers/left_primers_rc.fasta"
params.right_primers_rc = "$baseDir/resources/primers/right_primers_rc.fasta"
params.amplicon_middle_positions = "$baseDir/resources/amplicon_bed/amplicons_bed.txt"

params.force_sites = "$baseDir/force_sites.vcf.gz"

// cutadapt
params.quality_cutoff = 25
params.minimum_length = 60
params.maximum_length = 300

params.humans_index_dir = "$baseDir/resources/rtn_files/humans"
params.humans_index_base = "humans_modified.fa"

params.numts_index_dir = "$baseDir/resources/rtn_files/numts"
params.numts_index_base = "Calabrese_Dayama_Smart_Numts_modified.fa"


// params.humans = "$baseDir/resources/rtn_files/humans_modified.fa"
// params.humans_amb = "$baseDir/resources/rtn_files/humans_modified.fa.amb"
// params.humans_ann = "$baseDir/resources/rtn_files/humans_modified.fa.ann"
// params.humans_bwt = "$baseDir/resources/rtn_files/humans_modified.fa.bwt"
// params.humans_pac = "$baseDir/resources/rtn_files/humans_modified.fa.pac"
// params.humans_sa = "$baseDir/resources/rtn_files/humans_modified.fa.sa"

// params.numts = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa"
// params.numts_amb = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.amb"
// params.numts_ann = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.ann"
// params.numts_bwt = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.bwt"
// params.numts_pac = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.pac"
// params.numts_sa = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.sa"

params.mtdna_database = "$baseDir/HelixMTdb_20200327_short.vcf.gz"

params.fdstools_library = "$baseDir/resources/fdstools/mtNG_lib2_v211-flank.txt"

// rtn
params.mapQ = 30

// fdstools
params.minimum = 5
params.num_threads = 6 
params.min_reads_filt = 2
params.min_abs = 5 
params.min_pct_of_max = 0 
params.min_pct_of_sum = 3 
params.allele_min_abs = 5 
params.allele_min_pct_of_max = 0 
params.allele_min_pct_of_sum = 3
 

// mutect2
params.detection_limit = 0.08
params.baseQ = 30
params.callable_depth = 6
params.initial_tumor_lod = 0
params.tumor_lod_to_emit = 0
params.native_pair_hmm_threads = 1
params.max_reads_per_aglignment_start = 0

params.python_script_remove_scb = "$baseDir/resources/scripts/remove_soft_clipped_bases.py"
params.python_script2 = "$baseDir/resources/scripts/python_empop.py"
params.python_script3 = "$baseDir/resources/scripts/python_coverage.py"
params.python_script4 = "$baseDir/resources/scripts/process_fdstools_output.py"
params.python_script5 = "$baseDir/resources/scripts/process_mutect2_output.py"
params.python_script6 = "$baseDir/resources/scripts/merge_fdstools_mutect2.py"

    // rm -r "$baseDir/work"
    // rm -r "$baseDir/results"
    // rm .nextflow.*

log_text = """\
         m t D N A - N i m a G e n  P I P E L I N E    
         ==========================================
         mtDNA reference genome         : ${params.reference}
         reads                          : ${params.reads}
         
         MERGING (with FLASH)
         min-overlap                    : $params.min_overlap # The minimum required overlap length between two reads to provide a confident overlap (default: 10bp) 
         max-overlap                    : $params.max_overlap # Maximum overlap length expected in approximately 90% of read pairs. 
         max-mismatch-density           : $params.max_mismatch_density # Maximum allowed ratio between the number of mismatched base pairs and the overlap length 
         allow-outies                   : enabled (hard coded) # If a read pair can be combined in both "innie" and "outie" orientations, the better-fitting one will be chosen.

         TRIMMING (with CUTADAPT) 
         quality-cutoff                 : $params.quality_cutoff # Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal.
         minimum-length                 : $params.minimum_length # Discard reads shorter than LEN. Default: 0
         maximum-length                 : $params.maximum_length # Discard reads longer than LEN. Default: no limit
         discard-untrimmed              : enabled (hard coded) # Discard reads that do not contain an adapter.

         NUMTs REMOVAL (with RTN)
         mapQ                           : $params.mapQ # Used to filter out reads assigned as NUMTs by RTN 

         VARIANT CALLING (with FDSTOOLS)
         tssv: minimum                  : $params.minimum # report only sequences with this minimum number of reads (default: 2) 
         tssv: num-threads              : $params.num_threads # number of worker threads to use (default: 1)
         samplestats: min-reads-filt    : $params.min_reads_filt # the minimum number of reads (default: 1)
         vis: min-abs                   : $params.min_abs # only show sequences with this minimum number of reads (default: 5)
         vis: min-pct-of-max            : $params.min_pct_of_max # for sample: only show sequences with at least this percentage of the number of reads of the highest allele of a marker
         vis: min-pct-of-sum            : $params.min_pct_of_sum # only show sequences with at least this percentage of the total number of reads of a marker (default: 0.0)
         vis: allele-min-abs            : $params.allele_min_abs # the minimum number of reads (default: 30)
         vis: allele-min-pct-of-max     : $params.allele_min_pct_of_max # the minimum percentage of reads w.r.t. the highest allele of the marker (default: 2.0)
         vis: allele-min-pct-of-sum     : $params.allele_min_pct_of_sum # the minimum percentage of reads w.r.t. the markers total number of reads (default: 1.5)

         VARIANT CALLING (with MUTECT2)
         detection_limit                : $params.detection_limit #
         baseQ                          : $params.baseQ # 
         callable-depth                 : $params.callable_depth #
         initial-tumor-lod              : $params.initial_tumor_lod #
         tumor_lod_to_emit              : $params.tumor_lod_to_emit #
         native_pair_hmm_threads        : $params.native_pair_hmm_threads #
         max_reads_per_aglignment_start : $params.max_reads_per_aglignment_start #

         outdir                         : ${params.outdir}
         """

log.info(log_text)

process p00_pipeline_parameters{
    publishDir "$params.outdir", mode: params.publish_dir_mode

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
    path("${reference}.*")

    script:
    """
    bwa index $reference 
    
    samtools faidx $reference 
    samtools dict $reference
    """
}

process p02_map_raw_fastq_p01 {
    tag "p02: bwa mem on $sample_id"
    publishDir "$params.outdir/p02_mapped_w_scb_bam", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(reads)
    path reference
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}_R1.sam"), path("${sample_id}_R2.sam"), emit: p02_raw_sam_ch
    tuple path("${sample_id}_R1.bam"), path("${sample_id}_R2.bam"), path("${sample_id}_R1.bam.bai"), path("${sample_id}_R2.bam.bai")

    script:
    """
    mv ${reads[0]} tmp.fastq.gz
    cutadapt -a ${params.adapter} -o ${reads[0]}  tmp.fastq.gz 

    bwa mem ${reference} ${reads[0]} > ${sample_id}_R1.sam 
    bwa mem ${reference} ${reads[1]} > ${sample_id}_R2.sam 
    samtools view -bS ${sample_id}_R1.sam | samtools sort -o ${sample_id}_R1.bam
    samtools view -bS ${sample_id}_R2.sam | samtools sort -o ${sample_id}_R2.bam
    samtools index ${sample_id}_R1.bam
    samtools index ${sample_id}_R2.bam 
    """
}

process p03_filter_softclipped_sam_p02 {
    tag "p03: removing scb from $sample_id"
    publishDir "$params.outdir/p03_mapped_wo_scb_bam", mode: 'copy', pattern: '*sorted.bam*'

    input:
    tuple val(sample_id), path(sam_r1), path(sam_r2)
    path python_script_remove_scb

    output:
    tuple val(sample_id), path("${sample_id}_R1_wo_scb.bam"), path("${sample_id}_R2_wo_scb.bam"), emit: p03_bam_files_wo_scb_ch
    tuple val(sample_id), path("${sample_id}_R1_wo_scb_sorted.bam"), path("${sample_id}_R2_wo_scb_sorted.bam"), path("${sample_id}_R1_wo_scb_sorted.bam.bai"), path("${sample_id}_R2_wo_scb_sorted.bam.bai")

    script:
    """ 
    cat ${sam_r1} | python $python_script_remove_scb > ${sample_id}_R1_wo_scb.sam
    cat ${sam_r2} | python $python_script_remove_scb > ${sample_id}_R2_wo_scb.sam
    samtools view -Sb ${sample_id}_R1_wo_scb.sam > ${sample_id}_R1_wo_scb.bam
    samtools view -Sb ${sample_id}_R2_wo_scb.sam > ${sample_id}_R2_wo_scb.bam
    samtools sort -o ${sample_id}_R1_wo_scb_sorted.bam ${sample_id}_R1_wo_scb.bam
    samtools sort -o ${sample_id}_R2_wo_scb_sorted.bam ${sample_id}_R2_wo_scb.bam
    samtools index ${sample_id}_R1_wo_scb_sorted.bam
    samtools index ${sample_id}_R2_wo_scb_sorted.bam
    """
}

process p04_convert_bam_2_fastq_p03 {
    tag "p04: convert BAM_wo_scb to FASTQ on $sample_id"
    // publishDir "$params.outdir/p04_wo_scb_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_r1_wo_scb), path(bam_r2_wo_scb)
    
    output:
     tuple val(sample_id), path("${sample_id}_R1_wo_scb.fastq"), path("${sample_id}_R2_wo_scb.fastq")

    script:
    """
    samtools fastq ${bam_r1_wo_scb} > ${sample_id}_R1_wo_scb.fastq    
    samtools fastq ${bam_r2_wo_scb} > ${sample_id}_R2_wo_scb.fastq   
    """
}


process p05_merge_fastq_p04 {
    tag "p05: flash on $sample_id"
    // publishDir "$params.outdir/p05_wo_scb_merged_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_r1_wo_scb), path(fastq_r2_wo_scb)
    
    output:
    tuple val(sample_id), path("${sample_id}_wo_scb_merged.extendedFrags.fastq")

    script:
    """
    flash ${fastq_r1_wo_scb} ${fastq_r2_wo_scb} -m $params.min_overlap -M $params.max_overlap -x $params.max_mismatch_density -O -o ${sample_id}_wo_scb_merged
    """
}

process p06_trim_merged_fastq_p05 {
    tag "p06: cutadapt on $sample_id"
    // publishDir "$params.outdir/p06_wo_scb_merged_trimmed_fastq", mode: 'copy'

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

process p07_map_merged_fastq_p01_p05 {
    tag "p07: bwa mem on ${fastq_wo_scb_merged.baseName}"
    publishDir "$params.outdir/p07_mapped_wo_scb_merged_bam", mode: 'copy', pattern: '*.bam*'

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
    // samtools depth -a -b $amplicon_middle_positions ${sample_id}.bam > ${sample_id}_coverage.txt
// samtools view -h ${sample_id}_tmp.bam | awk '\$1 ~ /^@/ || \$6 !~ /S/' | samtools view -b -o ${sample_id}.bam

process p08_map_merged_trimmed_fastq_p01_p06 {
    tag "p08: bwa mem on ${fastq_wo_scb_merged_trimmed.baseName}"
    publishDir "$params.outdir/p08_mapped_wo_scb_merged_trimmed_bam", mode: 'copy', pattern: '*.bam*'

    input:
    tuple val(sample_id), path(fastq_wo_scb_merged_trimmed)
    path reference
    path index_files

    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${fastq_wo_scb_merged_trimmed.baseName}.bam"), path("${fastq_wo_scb_merged_trimmed.baseName}.bam.bai"), path("${fastq_wo_scb_merged_trimmed.baseName}_coverage.txt")

    script:
    """
    bwa mem ${reference} ${fastq_wo_scb_merged_trimmed} | samtools view -Sb - | samtools sort -o ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam

    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${fastq_wo_scb_merged_trimmed.baseName}_sorted.bam  -o ${fastq_wo_scb_merged_trimmed.baseName}.bam
    samtools index ${fastq_wo_scb_merged_trimmed.baseName}.bam

    samtools depth -a -b $amplicon_middle_positions ${fastq_wo_scb_merged_trimmed.baseName}.bam > ${fastq_wo_scb_merged_trimmed.baseName}_coverage.txt
    """
}

process p09_filter_numts_merged_bam_p07 {
    tag "p09: rtn on ${bam_wo_scb_merged.baseName}"
    publishDir "$params.outdir/p09_filtered_numts_bam_for_fdstoold", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_wo_scb_merged), path(bam_wo_scb_merged_index)
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
    samtools view -h -q $params.mapQ ${bam_wo_scb_merged.baseName}.rtn.bam > ${bam_wo_scb_merged.baseName}.rtn_tmp.bam
    samtools fastq ${bam_wo_scb_merged.baseName}.rtn_tmp.bam > ${bam_wo_scb_merged.baseName}_wo_NUMTs.fastq
    """
}
    // samtools view -h -q $params.mapQ ${sample_id}.rtn.bam > ${sample_id}.rtn_tmp.bam
    // samtools depth -a -b $amplicon_middle_positions ${sample_id}.rtn_tmp.bam > ${sample_id}_coverage_numts.txt

process p10_filter_numts_trimmed_merged_bam_p08 {
    tag "p10: rtn on ${bam_wo_scb_merged_trimmed.baseName}"
    publishDir "$params.outdir/p10_filtered_numts_bam_for_mutect2", mode: 'copy', pattern: '*.bam*'
    publishDir "$params.outdir/p10_coverage_txt", mode: 'copy', pattern: '*.txt'

    input:
    tuple val(sample_id), path(bam_wo_scb_merged_trimmed), path(bam_wo_scb_merged_trimmed_index), path(coverage_txt)
    path amplicon_middle_positions
    path humans_index
    val humans_base
    path numts_index
    val numts_base


    output:
    tuple val(sample_id), path("${bam_wo_scb_merged_trimmed.baseName}.rtn.bam"), path("${bam_wo_scb_merged_trimmed.baseName}.rtn.bam.bai"), path(coverage_txt), path("${bam_wo_scb_merged_trimmed.baseName}_coverage_wo_NUMTs.txt")

    script:
    """
    rtn -h "${humans_index}/${humans_base}" -n "${numts_index}/${numts_base}" -b $bam_wo_scb_merged_trimmed
    samtools view -h -q $params.mapQ ${bam_wo_scb_merged_trimmed.baseName}.rtn.bam > ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam
    samtools depth -a -b $amplicon_middle_positions ${bam_wo_scb_merged_trimmed.baseName}_filtered.rtn.bam > ${bam_wo_scb_merged_trimmed.baseName}_coverage_wo_NUMTs.txt
    """
}
    // samtools fastq ${sample_id}.rtn_tmp.bam > ${sample_id}_rtn_filtered.fastq
    // samtools fastq ${sample_id}.rtn_tmp.bam > ${sample_id}_rtn_filtered.fastq
    // rtn -p -h $humans -n $numts -b $trimmed_bam_file

process p11_variant_calling_fdstools_p09 {
    tag "ka: fdstools on $sample_id"
    publishDir "$params.outdir/ka_fdstools", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(rtn_fastq)

    output:
    tuple  val(sample_id), path("${sample_id}.tssv.csv"), path("${sample_id}.report.txt"), path("${sample_id}.sc.csv"), path("${sample_id}.sast.csv"), path("${sample_id}.html")

    """
    fdstools tssv $params.fdstools_library  ${rtn_fastq} ${sample_id}.tssv.csv --minimum $params.minimum --num-threads $params.num_threads --report ${sample_id}.report.txt
	fdstools seqconvert allelename ${sample_id}.tssv.csv ${sample_id}.sc.csv --library $params.fdstools_library
	fdstools samplestats --min-reads-filt $params.min_reads_filt ${sample_id}.sc.csv ${sample_id}.sast.csv
	fdstools vis --min-abs $params.min_abs --min-pct-of-max $params.min_pct_of_max --min-pct-of-sum $params.min_pct_of_sum --allele-min-abs $params.allele_min_abs --allele-min-pct-of-max $params.allele_min_pct_of_max --allele-min-pct-of-sum $params.allele_min_pct_of_sum sample ${sample_id}.sast.csv ${sample_id}.html
    """
}


process i_CALCULATE_STATISTICS {
    tag "i: calculate_statistics on $sample_id"
    publishDir "$params.outdir/i_calculate_statistics", mode: 'copy'

    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt), path(coverage_numts_txt), path(rtn_fastq)
    path python_script3

    output:
    path "*summary.txt", emit: stats_ch
    path "*mapping.txt", emit: mapping_ch
    path "*.zip", emit: fastqc_ch
    path("*.bam"), includeInputs: true, emit: fixed_file
    path("*coverage_plot.png")

    script:
    def output_name = "${sample_id}.summary.txt"
    def mapping_name = "${sample_id}.mapping.txt"
 
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    
    // 16623 - lenght of the extended rCRS
    // samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${bam_file}  -o ${sample_id}_tmp.bam
    // mv ${sample_id}_tmp.bam ${bam_file}

    """
    ## Create Mapping File
    ## mkdir "${sample_id}"
    ## echo -e "Sample\tFilename" > /${sample_id}/$mapping_name
    echo -e "Sample\tFilename" > $mapping_name

    echo "\$(samtools samples ${bam_file})" >> $mapping_name

    ## Calculate summary statistics
    samtools coverage ${bam_file} > samtools_coverage_${sample_id}.txt
    csvtk grep -t -f3 -p 16623 -C '\$' samtools_coverage_${sample_id}.txt -T -o mtdna.txt --num-cpus ${task.cpus} 
        
    contig=\$(csvtk cut -t -f 1 mtdna.txt --num-cpus ${task.cpus})
    numreads=\$(csvtk cut -t -f 4 mtdna.txt --num-cpus ${task.cpus})
    covered_bases=\$(csvtk cut -t -f 5 mtdna.txt --num-cpus ${task.cpus})
    covered_bases_percentage=\$(csvtk cut -t -f 6 mtdna.txt --num-cpus ${task.cpus})
    mean_depth=\$(csvtk cut -t -f 7 mtdna.txt --num-cpus ${task.cpus})
    mean_base_quality=\$(csvtk cut -t -f 8 mtdna.txt --num-cpus ${task.cpus})
    mean_map_quality=\$(csvtk cut -t -f 9 mtdna.txt --num-cpus ${task.cpus})
    readgroup=\$(samtools view -H ${bam_file} | csvtk grep -I -H -r -p "^@RG" --num-cpus ${task.cpus} | sed 's/\t/,/g' | head -n 1)
    
    echo -e "Sample\tParameter\tValue" > $output_name
    echo -e "${bam_file}\tContig\t\${contig}" >> $output_name
    echo -e "${bam_file}\tNumberofReads\t\${numreads}" >> $output_name
    echo -e "${bam_file}\tCoveredBases\t\${covered_bases}" >> $output_name
    echo -e "${bam_file}\tCoveragePercentage\t\${covered_bases_percentage}" >> $output_name
    echo -e "${bam_file}\tMeanDepth\t\${mean_depth}" >> $output_name
    echo -e "${bam_file}\tMeanBaseQuality\t\${mean_base_quality}" >> $output_name
    echo -e "${bam_file}\tMeanMapQuality\t\${mean_map_quality}" >> $output_name
    echo -e "${bam_file}\tRG\t\${readgroup}" >> $output_name

    fastqc --threads ${task.cpus} --memory ${avail_mem} $bam_file -o .
    
    sleep 2
    rm -f ${sample_id}_coverage_plot.png
    python $python_script3 $coverage_txt $coverage_numts_txt ${sample_id}_coverage_plot.png

    """
}
// echo "\$(${sample_id} ${bam_file})" >> $mapping_name

process j_INDEX_CREATION {
	tag "j: samtools on $reference"
    publishDir "$params.outdir/j_index_creation", mode: 'copy'
    
    input:
	path reference
	val mtdna_tag

	output:
	path "ref*.{dict,fai}", emit: fasta_index_ch
	path "ref.fasta", emit: ref_ch

	"""
	sed -e "s/^>.*/>$mtdna_tag/" $reference > ref.fasta
    samtools faidx ref.fasta
    	samtools dict ref.fasta \
	    -o ref.dict
	"""
}

process p12_variant_calling_mutect2_p10 {
    tag "p12: mutect2 on $sample_id"
    publishDir "$params.outdir/kb_mutect2", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt), path(coverage_txt_numts), path(rtn_fastq)
    path reference
    path fasta_index_files
    val detected_contig
    val method

    output:
    tuple  val(sample_id), path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), val(method), emit: mutect2_ch
    // path("${bam_file.baseName}.bamout")
    path("${bam_file.baseName}_sorted.bamout.bam")
    path("${bam_file.baseName}_sorted.bamout.bam.bai")
    // path("run1.vcf.gz ")
    // path("run2.vcf.gz ")
    // path("${sample_id}_mpileup.vcf.gz")
    // path("${sample_id}_mpileup_bamout.vcf.gz")

    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    
    // samtools index ${bam_file}
    
    """  
    mkdir tmp_${sample_id}


    gatk --java-options "-Xmx16G" \
        Mutect2 \
        -R ${reference} \
        -L 'chrM' \
        --min-base-quality-score ${params.baseQ} \
        --callable-depth $params.callable_depth \
        --linked-de-bruijn-graph true \
        --recover-all-dangling-branches true \
        --initial-tumor-lod 0   \
        --tumor-lod-to-emit 0   \
        --native-pair-hmm-threads 1 \
        --max-reads-per-alignment-start 0 \
        --bam-output ${bam_file.baseName}.bamout.bam \
        --tmp-dir tmp_${sample_id} \
        -I ${bam_file} \
        -O raw.vcf.gz 



    samtools sort -o ${bam_file.baseName}_sorted.bamout.bam ${bam_file.baseName}.bamout.bam
    samtools index ${bam_file.baseName}_sorted.bamout.bam 

    gatk  --java-options "-Xmx16G" \
        FilterMutectCalls \
        -R ${reference} \
        --min-reads-per-strand 0 \
        -V raw.vcf.gz \
        --tmp-dir . \
        -O ${bam_file.baseName}.vcf.gz

    bcftools norm \
        -m-any \
        -f ${reference} \
        -o ${bam_file.baseName}.norm.vcf.gz -Oz \
        ${bam_file.baseName}.vcf.gz 

    bcftools view \
    -i 'FORMAT/AF>=${params.detection_limit}' \
    -o ${bam_file.baseName}.vcf.gz -Oz \
    ${bam_file.baseName}.norm.vcf.gz 
    
    tabix -f ${bam_file.baseName}.vcf.gz

    rm ${bam_file.baseName}.norm.vcf.gz 
    
    """
}
    //     gatk --java-options "-Xmx16G" \
    //     Mutect2 \
    //     -R ${reference} \
    //     -L '${detected_contig}' \
    //     --min-base-quality-score 30 \
    //     --callable-depth 6 \
    //     --initial-tumor-lod -10   \
    //     --tumor-lod-to-emit -10   \
    //     --native-pair-hmm-threads 1 \
    //     --max-reads-per-alignment-start 0 \
    //     --bam-output ${bam_file.baseName}_1.bamout.bam \
    //     --tmp-dir tmp_${sample_id} \
    //     -I ${bam_file} \
    //     -O run1.vcf.gz 


    // bcftools isec -C -w1 -Oz -o run2_unique.vcf.gz run2.vcf.gz run1.vcf.gz
    // tabix -p vcf run2_unique.vcf.gz 
    // bcftools concat -a -Oz -o merged.vcf.gz run1.vcf.gz run2_unique.vcf.gz
    // bcftools sort merged.vcf.gz -Oz -o raw.vcf.gz
    // tabix -f raw.vcf.gz
    // cp run1.vcf.gz.stats raw.vcf.gz.stats


    // samtools sort -o ${bam_file.baseName}_sorted.bamout.bam ${bam_file.baseName}_1.bamout.bam
    // samtools index ${bam_file.baseName}_sorted.bamout.bam 


        // --genotype-pon-sites true \
        // --panel-of-normals $params.force_sites \

// --mitochondria-mode true \
// gatk --java-options "-Xmx16G"    \
//       Mutect2  \
//       -R ref.fasta   \
//       -L 'chrM'   \
//       --min-base-quality-score 32  \
//       --callable-depth 2  \
//       --initial-tumor-lod -10   \
//       --tumor-lod-to-emit -10   \
//       --linked-de-bruijn-graph true         --max-reads-$

// --bam-output s24-12889a_S3_L001.rtn.bamout.bam         
// --tmp-dir .         
// -I s24-12889a_S3_L001.rtn.bam         
// -O raw.vcf.gz
        //         --mitochondria-mode true \

        // --genotype-germline-sites true \
        // --linked-de-bruijn-graph true \
        // --mitochondria-mode true \

        // --af-of-alleles-not-in-resource 1e-3 \
        // --germline-resource ${params.mtdna_database} \

//  samtools view -H ${bam_file.baseName}.bamout.bam | grep -v "HaplotypeBAMWriter" > clean_header.sam
//     samtools reheader clean_header.sam ${bam_file.baseName}.bamout.bam > ${bam_file.baseName}.bamout.clean.bam

//     samtools sort -o ${bam_file.baseName}_sorted.bamout.clean.bam ${bam_file.baseName}.bamout.clean.bam
//     samtools index ${bam_file.baseName}_sorted.bamout.clean.bam

//     gatk --java-options "-Xmx16G" \
//         Mutect2 \
//         -R ${reference} \
//         -L '${detected_contig}' \
//         --min-base-quality-score ${params.baseQ} \
//         --callable-depth 2 \
//         --initial-tumor-lod 0 \
//         --tumor-lod-to-emit 0 \
//         --bam-output ${bam_file.baseName}.bamout2.bam \
//         --genotype-germline-sites true \
//         --linked-de-bruijn-graph true \
//         --max-reads-per-alignment-start 0 \
//         --af-of-alleles-not-in-resource 1e-3 \
//         --germline-resource ${params.mtdna_database} \
//         --tmp-dir . \
//         -I ${bam_file.baseName}_sorted.bamout.clean.bam \
//         -O ${bam_file.baseName}.vcf.gz



//     bcftools mpileup \
//         -f ${reference} \
//         --redo-BAQ \
//         --min-MQ 30 \
//         --min-BQ 32 \
//         -d 4000 \
//         -L 4000 \
//         -Ou \
//         ${bam_file} | \
//     bcftools call \
//         --ploidy 1 \
//         --multiallelic-caller \
//         --variants-only \
//         -Oz -o ${sample_id}_mpileup.vcf.gz






        // --max-reads-per-alignment-start 0 \

    
    // bcftools mpileup \
    //     -f ${reference} \
    //     --min-BQ 30 \
    //     -d 4000 \
    //     -L 4000 \
    //     -Ou \
    //     ${bam_file.baseName}.bamout.bam | \
    // bcftools call \
    //     --ploidy 1 \
    //     --multiallelic-caller \
    //     --variants-only \
    //     -Oz -o ${sample_id}_mpileup_bamout.vcf.gz

    // bcftools mpileup \
    //     --redo-BAQ \
    //     --min-BQ ${params.baseQ} \
    //     --per-sample-mF \
    //     --ploidy 1 \
    //     -f "${reference}" \
    //     --BCF ${bam_file} | \
    // bcftools call \
    //     --multiallelic-caller \
    //     --variants-only \
    //     -Oz > ${bam_file.baseName}_mpileup.vcf.gz


    // bcftools mpileup \
    //     --redo-BAQ \
    //     --min-BQ ${params.baseQ} \
    //     --per-sample-mF \
    //     -f "${reference}" \
    //     --BCF ${bam_file} | \
    // bcftools call \
    //     --multiallelic-caller \
    //     --variants-only \
    //     -Oz > ${bam_file.baseName}_mpileup.vcf.gz

    //  --output-tags DP,AD \
// --linked-de-bruijn-graph true \ TEST THIS TOMORROW! 27. March 2025
    // gatk  --java-options "-Xmx16G" \
    //     FilterMutectCalls \
    //     -R ${reference} \
    //     --min-reads-per-strand 0 \
    //     -V raw.vcf.gz \
    //     --tmp-dir . \
    //     -O ${bam_file.baseName}.vcf.gz

        // --genotype-germline-sites true \
        // --dont-use-soft-clipped-base true \
        // --linked-de-bruijn-graph true \
    // --normal-lod -20 \
        // --mitochondria-mode true \
        // --native-pair-hmm-threads 7 \
        // --max-reads-per-alignment-start 0 \
        // --pruning-lod-threshold -5 \
        // --af-of-alleles-not-in-resource 0.01 \
// --mitochondria-mode true \
// --af-of-alleles-not-in-resource 0.000001 \
// rm raw.vcf.gz

    //   --mitochondria-mode \
    //     --initial-tumor-lod 0 \
    //     --tumor-lod-to-emit 0 \


// /Users/peter/anaconda3/pkgs/gatk4-4.6.1.0-py310hdfd78af_0/share/gatk4-4.6.1.0-0/gatk  --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \
// /Users/peter/anaconda3/pkgs/gatk4-4.6.1.0-py310hdfd78af_0/share/gatk4-4.6.1.0-0/gatk  --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \

process l_FINAL_VARIANTS {
    tag "l: final_variants on $sample_id"
    publishDir "$params.outdir/l_final_variants", mode: 'copy'
    // publishDir "${params.output}", mode: 'copy'

    input:
    tuple  val(sample_id), path(vcf_file), path(vcf_file_idx), val(method), path(tssv_file), path(report_file), path(sc_file), path(sast_file), path(html_file)
    // tuple  val(sample_id), , 
    path reference
    path python_script2
    path python_script4
    path python_script5

    output:
    tuple val(sample_id), path("${vcf_file.baseName}.${method}.filtered.empop.txt"), path("${vcf_file.baseName}.${method}.filtered.empop_final.txt"), path("${sample_id}_fdstools_processed.txt"), emit: fdstools_mutect2_variants_ch

    script:
    def vcf_name = "${vcf_file}".replaceAll('.vcf.gz', '')

    """
    echo -e "test"
    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tMeanBaseQuality\tCoverage\tGT" \
        > ${vcf_file.baseName}.${method}.txt

    bcftools query -u \
        -f "${vcf_name}.bam\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%MBQ\t%AD\t%GT]\n" \
        ${vcf_file} | \
        awk -F'\t' '(\$2 !~ /bla/)' \
        >> ${vcf_file.baseName}.${method}.txt
         
    

    ## annotating SNVS and INDELs for reporting
    awk 'BEGIN {OFS="\t"} {
        if (NR == 1) { print \$0, "Type"; next }
        if ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5)) { \$10="3" }
        else { \$10="2" }
        print
    }' ${vcf_file.baseName}.${method}.txt > ${vcf_file.baseName}.${method}.filtered.txt

    python $python_script2 ${vcf_file.baseName}.${method}.filtered.txt ${vcf_file.baseName}.${method}.filtered.empop.txt $reference 
    python $python_script5 ${vcf_file.baseName}.${method}.filtered.empop.txt ${vcf_file.baseName}.${method}.filtered.empop_final.txt
    python $python_script4 ${sast_file} ${sample_id}_fdstools_processed.txt 
    """
}
        // else if (\$9 == "0/1" || \$9 == "1/0" || \$9 == "0|1" || \$9 == "1|0") { \$10="2" }


// process m_PROCESS_FDSTOOLS {
//     tag "m: process_fdstools on $sample_id"
//     publishDir "$params.outdir/m_processed_fdstools", mode: 'copy'
//     // publishDir "${params.output}", mode: 'copy'

//     input:
//     tuple  val(sample_id), path(tssv_file), path(report_file), path(sc_file), path(sast_file), path(html_file)

//     output:
//     tuple val(sample_id), path("${sample_id}_fdstools_processed.txt"), emit: fdstools_variants_ch


//     """
//     python $python_script4 ${sast_file} ${sample_id}_fdstools_processed.txt 
//     """
// }

process n_MERGE_FDSTOOLS_MUTECT2 {
    tag "n: merged_callers on $sample_id"
    publishDir "$params.outdir/n_merged_callers", mode: 'copy'
    // publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_mutect2), path(filtered_mutect2_final), path(filtered_fdstools)
    path python_script6
    // tuple val(sample_id), path(filtered_fdstools)

    output:
    path("${sample_id}_merged_variants.txt"), emit: merged_variants_ch


    """
    python $python_script6 ${filtered_fdstools} ${filtered_mutect2_final} ${sample_id}_merged_variants.txt

    """
}




workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    p00_pipeline_parameters(log_text)
 
    index_ch = p01_index_reference_fasta(params.reference)
  
    def detected_contig = "chrM"

    p02_map_raw_fastq_p01(read_pairs_ch, params.reference, index_ch)
    mapping_ch = p02_map_raw_fastq_p01.out.p02_raw_sam_ch
    
    p03_filter_softclipped_sam_p02(mapping_ch, params.python_script_remove_scb)
    cleaned_ch = p03_filter_softclipped_sam_p02.out.p03_bam_files_wo_scb_ch
    // cleaned_ch.waitFor()

    fastq_ch = p04_convert_bam_2_fastq_p03(cleaned_ch)
    // fastq_ch.waitFor()

    merging_ch = p05_merge_fastq_p04(fastq_ch)
    // merging_ch.waitFor()

    trimming_ch = p06_trim_merged_fastq_p05(merging_ch, params.left_primers, params.right_primers_rc)
    // trimming_ch.waitFor()

    mapping_final_ch = p07_map_merged_fastq_p01_p05(merging_ch, params.reference, index_ch, params.amplicon_middle_positions)
    // // mapping_final_ch.waitFor()

    mapping2_final_ch = p08_map_merged_trimmed_fastq_p01_p06(trimming_ch, params.reference, index_ch, params.amplicon_middle_positions)
    // mapping_final_ch.waitFor()

    humans_index_ch = Channel.value(file(params.humans_index_dir))
    numts_index_ch  = Channel.value(file(params.numts_index_dir))

    humans_base_ch = Channel.value(params.humans_index_base)
    numts_base_ch  = Channel.value(params.numts_index_base)

    rtn_ch = p09_filter_numts_merged_bam_p07(mapping_final_ch, params.amplicon_middle_positions, humans_index_ch, humans_base_ch, numts_index_ch, numts_base_ch)
    // // rtn_ch.waitFor()

    rtn2_ch = p10_filter_numts_trimmed_merged_bam_p08(mapping2_final_ch, params.amplicon_middle_positions,  humans_index_ch, humans_base_ch, numts_index_ch, numts_base_ch)
    // // rtn_ch.waitFor()


    // i_ch = i_CALCULATE_STATISTICS(rtn_ch, params.python_script3)
    // // i_ch.waitFor()

    // j_ch = j_INDEX_CREATION(params.reference, detected_contig )
    // // j_ch.waitFor()

    // k_ch = k_MUTECT2(rtn_ch, j_INDEX_CREATION.out.ref_ch, j_INDEX_CREATION.out.fasta_index_ch, detected_contig, "mutect2_fusion")
    // // k_ch.waitFor()

    // vcf_ch = k_MUTECT2.out.mutect2_ch 

    // // l_FINAL_VARIANTS(vcf_ch, params.reference, params.python_script2)
    
    fdstools_ch = p11_variant_calling_fdstools_p09(rtn_ch)

    // final_inputs = vcf_ch.join(fdstools_ch, by: 0)
    // //     .map { sid, vcf_parts, fdstools_parts -> tuple(sid, *vcf_parts, *fdstools_parts) }

    // // final_inputs = vcf_ch
    // // .join(fdstools_ch, by: 0)
    // // .map { sid, vcf_parts, fdstools_parts ->
    // //     tuple(sid, *vcf_parts, *fdstools_parts)  // <-- this flattens the list properly
    // // }

    // // final_inputs.view()
    // // final_inputs.into { final_variants_inputs }

    // l_ch = l_FINAL_VARIANTS(final_inputs, params.reference, params.python_script2, params.python_script4, params.python_script5)
    // variant_callers_ch = l_FINAL_VARIANTS.out.fdstools_mutect2_variants_ch

    // n_ch = n_MERGE_FDSTOOLS_MUTECT2(variant_callers_ch, params.python_script6)

}