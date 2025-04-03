nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/raw_data/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/resources/rCRS/rCRS2.fasta"
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 140 // default in FLASH is 65
params.max_mismatch_density = 0.25 //default in FLASH is 0.25
// params.multiqc = "$baseDir/multiqc"
params.publish_dir_mode = "symlink"
params.outdir = "results.nosync"

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

params.humans = "$baseDir/resources/rtn_files/humans.fa"
params.humans_amb = "$baseDir/resources/rtn_files/humans.fa.amb"
params.humans_ann = "$baseDir/resources/rtn_files/humans.fa.ann"
params.humans_bwt = "$baseDir/resources/rtn_files/humans.fa.bwt"
params.humans_pac = "$baseDir/resources/rtn_files/humans.fa.pac"
params.humans_sa = "$baseDir/resources/rtn_files/humans.fa.sa"

params.numts = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa"
params.numts_amb = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.amb"
params.numts_ann = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.ann"
params.numts_bwt = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.bwt"
params.numts_pac = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.pac"
params.numts_sa = "$baseDir/resources/rtn_files/Calabrese_Dayama_Smart_Numts_modified.fa.sa"

params.mtdna_database = "$baseDir/HelixMTdb_20200327_short.vcf.gz"

params.fdstools_library = "$baseDir/resources/fdstools/mtNG_lib2_v211-flank.txt"


// mutect2
params.detection_limit = 0.08
params.mapQ = 30
params.baseQ = 32
params.alignQ = 30

params.python_script = "$baseDir/resources/scripts/remove_soft_clipped_bases.py"
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
         mtDNA reference genome   : ${params.reference}
         reads                    : ${params.reads}
         
         MERGING (with FLASH)
         merging_min-overlap      : $params.min_overlap
         merging_max-overlap      : $params.max_overlap
         max-mismatch-density     : $params.max_mismatch_density
         
         TRIMMING (with CUTADAPT) 
         quality_cutoff           : $params.quality_cutoff
         minimum_length           : $params.minimum_length
         maximum_length           : $params.maximum_length

         VARIANT CALLING (with MUTECT2)
         detection_limit          : $params.detection_limit
         mapQ                     : $params.mapQ
         baseQ                    : $params.baseQ
         alignQ                   : $params.alignQ

         outdir                   : ${params.outdir}
         """

//  publish_dir_mode       : $params.publish_dir_mode
// read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

// read_pairs_ch.view()

log.info(log_text)

// Save parameters into file
process a_write_log{
    publishDir "$params.outdir", mode: params.publish_dir_mode

    input:
    val(logs)

    output:
    path "parameters.txt"

    script:
    """
    echo '$logs' > parameters.txt
    """
}

process b_INDEXING {
    tag "b: bwa index on $reference"   
    // publishDir "$params.outdir/b_index", mode: 'copy'

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
//     // sed -e "s/^>.*/>$mtdna_tag/" $reference > ref.fasta

process c_MAPPING_2_SAM {
    tag "c: bwa mem on $sample_id"
    publishDir "$params.outdir/c_mapped_sam", mode: 'copy'

    input:
    path reference
    path index_files
    tuple val(sample_id), path(reads)
    // path merged_file
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.sam"), path("${sample_id}_R2.sam"), emit: mapping_test
    // , emit: mapping_sam_ch
    tuple path("${sample_id}_R1.bam"), path("${sample_id}_R2.bam"), path("${sample_id}_R1.bam.bai"), path("${sample_id}_R2.bam.bai"), emit: test

    script:
    """
    mv ${reads[0]} tmp.fastq.gz
    cutadapt -a ATCATAACAAAAAATTTCCACCAAA -o ${reads[0]}  tmp.fastq.gz 

    bwa mem $reference ${reads[0]} > ${sample_id}_R1.sam 
    bwa mem $reference ${reads[1]} > ${sample_id}_R2.sam 
    samtools view -bS ${sample_id}_R1.sam | samtools sort -o ${sample_id}_R1.bam
    samtools view -bS ${sample_id}_R2.sam | samtools sort -o ${sample_id}_R2.bam
    samtools index ${sample_id}_R1.bam
    samtools index ${sample_id}_R2.bam 
    """
}
    // cutadapt -a ATCATAACAAAAAATTTCCACCAAA -o tmp2_R1.fastq.gz  ${reads[0]} 

    // cutadapt -a file:$right_primers_rc  -o tmp_R1.fastq.gz tmp2_R1.fastq.gz -j 4
	// cutadapt -a file:$left_primers_rc -o tmp_R2.fastq.gz ${reads[1]}  -j 4


    // cutadapt -a file:$right_primers_rc  -o tmp_R1.fastq.gz tmp2_R1.fastq.gz -j 4
	// cutadapt -a file:$left_primers_rc -o tmp_R2.fastq.gz ${reads[1]}  -j 4
	
    // mv ${reads[0]} tmp.fastq.gz
    // cutadapt -a ATCATAACAAAAAATTTCCACCAAA -o ${reads[0]}  tmp.fastq.gz 

    // mv ${reads[0]} tmp.fastq.gz
    // cutadapt -a $params.adapter -o ${reads[0]}  tmp.fastq.gz 

process d_REMOVE_SOFT_CLIPPED_BASES {
    // container 'peterresutik/nimagen-pipeline:latest'
    tag "d: remove_scb on $sample_id"
    publishDir "$params.outdir/d_cleaned", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_r1), path(sam_r2)
    path script

    output:
    tuple val(sample_id), path("${sample_id}_R1_cleaned.bam"), path("${sample_id}_R2_cleaned.bam")
    // , path("${sample_id}_R1_cleaned.bam.bai"), path("${sample_id}_R2_cleaned.bam.bai")

    script:
    """ 
    cat ${sam_r1} | python $script > ${sample_id}_R1_cleaned.sam
    cat ${sam_r2} | python $script > ${sample_id}_R2_cleaned.sam
    samtools view -Sb ${sample_id}_R1_cleaned.sam > ${sample_id}_R1_cleaned.bam
    samtools view -Sb ${sample_id}_R2_cleaned.sam > ${sample_id}_R2_cleaned.bam

    """
}

    // samtools view -Sb ${sample_id}_R1_cleaned.sam | samtools sort -o ${sample_id}_R1_cleaned.bam
    // samtools view -Sb ${sample_id}_R2_cleaned.sam | samtools sort -o ${sample_id}_R2_cleaned.bam
    
    // samtools index ${sample_id}_R1_cleaned.bam
    // samtools index ${sample_id}_R2_cleaned.bam

process e_BACK_2_FASTQ {
    tag "e: convert_2_fastq on $sample_id"
    // publishDir "$params.outdir/e_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(cleaned_bam_r1), path(cleaned_bam_r2)
    // , path(cleaned_bai_r1), path(cleaned_bai_r2)
    // path merged_file
    
    output:
     tuple val(sample_id), path("${sample_id}_R1_cleaned.fastq"), path("${sample_id}_R2_cleaned.fastq")

    script:
    """
    samtools fastq ${cleaned_bam_r1} > ${sample_id}_R1_cleaned.fastq    
    samtools fastq ${cleaned_bam_r2} > ${sample_id}_R2_cleaned.fastq   
    """
}


process f_MERGING {
    tag "f: flash on $sample_id"
    // publishDir "$params.outdir/f_merged", mode: 'copy'

    input:
    tuple val(sample_id), path(cleaned_fastq_r1), path(cleaned_fastq_r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_cleaned_merged.extendedFrags.fastq")

    script:
    """
    flash ${cleaned_fastq_r1} ${cleaned_fastq_r2} -m $params.min_overlap -M $params.max_overlap -x $params.max_mismatch_density -O -o ${sample_id}_cleaned_merged
    """
}

process g_TRIMMING {
    tag "g: cutadapt on $sample_id"
    // publishDir "$params.outdir/g_cutadapt", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fastq)
    path left_primers
    path right_primers

    output:
    tuple val(sample_id), path("${sample_id}_cleaned_merged_trimmed_left_right.fastq")

    script:
    """    
    cutadapt -g file:$left_primers -q $params.quality_cutoff -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_cleaned_merged_trimmed_left.fastq $merged_fastq  
    cutadapt -a file:$right_primers -q $params.quality_cutoff -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_cleaned_merged_trimmed_left_right.fastq ${sample_id}_cleaned_merged_trimmed_left.fastq    
    """
}
// --discard-untrimmed

process ha_MAPPING_2_BAM_trimmed {
    tag "ha: bwa mem on $sample_id"
    publishDir "$params.outdir/ha_mapped_final", mode: 'copy'

    input:
    path reference
    path index_files
    // tuple val(sample_id), path(merged_fastq)
    tuple val(sample_id), path(trimmed_fastq)
    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}_coverage.txt")

    script:
    """
    bwa mem $reference ${trimmed_fastq} | samtools view -Sb - > ${sample_id}_tmp.bam    
    samtools view -h ${sample_id}_tmp.bam | awk '\$1 ~ /^@/ || \$6 !~ /S/' | samtools view -b -o ${sample_id}.bam
    samtools sort -o ${sample_id}_tmp.bam ${sample_id}.bam
    
    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${sample_id}_tmp.bam  -o ${sample_id}_tmp2.bam
    mv ${sample_id}_tmp2.bam ${sample_id}.bam
    samtools index ${sample_id}.bam

    samtools depth -a -b $amplicon_middle_positions ${sample_id}.bam > ${sample_id}_coverage.txt

    """
}
    // bwa mem $reference ${trimmed_fastq} | samtools view -Sb - > ${sample_id}_trimmed_tmp3.bam    
    // samtools view -h ${sample_id}_trimmed_tmp3.bam | awk '\$1 ~ /^@/ || \$6 !~ /S/' | samtools view -b -o ${sample_id}_trimmed_tmp4.bam
    // samtools sort -o ${sample_id}_trimmed_tmp.bam ${sample_id}_trimmed_tmp4.bam
    
    // samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${sample_id}_trimmed_tmp.bam  -o ${sample_id}_trimmed_tmp2.bam
    // mv ${sample_id}_trimmed_tmp2.bam ${sample_id}_trimmed.bam
    // samtools index ${sample_id}_trimmed.bam


process ha2_MAPPING_2_BAM_merged {
    tag "ha2: bwa mem on $sample_id"
    publishDir "$params.outdir/ha_mapped_final", mode: 'copy'

    input:
    path reference
    path index_files
    tuple val(sample_id), path(merged_fastq)
    // tuple val(sample_id), path(trimmed_fastq)
    path amplicon_middle_positions

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}_coverage.txt")

    script:
    """
    bwa mem $reference ${merged_fastq} | samtools view -Sb - > ${sample_id}_tmp.bam    
    samtools view -h ${sample_id}_tmp.bam | awk '\$1 ~ /^@/ || \$6 !~ /S/' | samtools view -b -o ${sample_id}.bam
    samtools sort -o ${sample_id}_tmp.bam ${sample_id}.bam
    
    samtools addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${sample_id}_tmp.bam  -o ${sample_id}_tmp2.bam
    mv ${sample_id}_tmp2.bam ${sample_id}.bam
    samtools index ${sample_id}.bam

    samtools depth -a -b $amplicon_middle_positions ${sample_id}.bam > ${sample_id}_coverage.txt

    """
}

process hb_NUMTs_trimmed {
    tag "hb: rtn on $sample_id"
    publishDir "$params.outdir/hb_numts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt)
    path amplicon_middle_positions
    path humans
    path humans_amb
    path humans_ann
    path humans_bwt
    path humans_pac
    path humans_sa
    path numts
    path numts_amb
    path numts_ann
    path numts_bwt
    path numts_pac
    path numts_sa

    output:
    tuple val(sample_id), path("${sample_id}.rtn.bam"), path("${sample_id}.rtn.bam.bai"), path(coverage_txt), path("${sample_id}_coverage_numts.txt"), path("${sample_id}_rtn_filtered.fastq")

    script:
    """
    rtn -p -h $humans -n $numts -b $bam_file

    samtools view -h -q 30 ${sample_id}.rtn.bam > ${sample_id}.rtn_tmp.bam

    samtools depth -a -b $amplicon_middle_positions ${sample_id}.rtn_tmp.bam > ${sample_id}_coverage_numts.txt
    
    samtools fastq ${sample_id}.rtn_tmp.bam > ${sample_id}_rtn_filtered.fastq

    """
}

process hb2_NUMTs_merged {
    tag "hb2: rtn on $sample_id"
    publishDir "$params.outdir/hb2_numts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt)
    path amplicon_middle_positions
    path humans
    path humans_amb
    path humans_ann
    path humans_bwt
    path humans_pac
    path humans_sa
    path numts
    path numts_amb
    path numts_ann
    path numts_bwt
    path numts_pac
    path numts_sa

    output:
    tuple val(sample_id), path("${sample_id}.rtn.bam"), path("${sample_id}.rtn.bam.bai"), path(coverage_txt), path("${sample_id}_coverage_numts.txt"), path("${sample_id}_rtn_filtered.fastq")

    script:
    """
    rtn -p -h $humans -n $numts -b $bam_file

    samtools view -h -q 30 ${sample_id}.rtn.bam > ${sample_id}.rtn_tmp.bam

    samtools depth -a -b $amplicon_middle_positions ${sample_id}.rtn_tmp.bam > ${sample_id}_coverage_numts.txt
    
    samtools fastq ${sample_id}.rtn_tmp.bam > ${sample_id}_rtn_filtered.fastq

    """
}

    // rtn -p -h $humans -n $numts -b $trimmed_bam_file

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

process FDSTOOLS{
    tag "ka: fdstools on $sample_id"
    publishDir "$params.outdir/ka_fdstools", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt), path(coverage_txt_numts), path(rtn_fastq)
    // path reference
    // path fasta_index_files
    // val detected_contig
    // val method

    output:
    tuple  val(sample_id), path("${sample_id}.tssv.csv"), path("${sample_id}.report.txt"), path("${sample_id}.sc.csv"), path("${sample_id}.sast.csv"), path("${sample_id}.html")

    """

    fdstools tssv $params.fdstools_library  ${rtn_fastq} ${sample_id}.tssv.csv --minimum 5 --num-threads 6 --report ${sample_id}.report.txt

	fdstools seqconvert allelename ${sample_id}.tssv.csv ${sample_id}.sc.csv --library $params.fdstools_library
	
	fdstools samplestats --min-reads-filt 2 ${sample_id}.sc.csv ${sample_id}.sast.csv
	
	fdstools vis --min-abs 5 --min-pct-of-max 0 --min-pct-of-sum 3 --allele-min-abs 5 --allele-min-pct-of-max 0 --allele-min-pct-of-sum 3 sample ${sample_id}.sast.csv ${sample_id}.html
    """
}


process k_MUTECT2 {
    tag "kb: mutect2 on $sample_id"
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
    gatk --java-options "-Xmx16G -XX:ActiveProcessorCount=1" \
        Mutect2 \
        -R ${reference} \
        -L '${detected_contig}' \
        --min-base-quality-score ${params.baseQ} \
        --callable-depth 2 \
        --linked-de-bruijn-graph true \
        --dont-use-soft-clipped-bases true \
        --alleles $params.force_sites \
        --initial-tumor-lod -10   \
        --tumor-lod-to-emit -10   \
        --max-reads-per-alignment-start 0 \
        --bam-output ${bam_file.baseName}.bamout.bam \
        --tmp-dir tmp_${sample_id} \
        -I ${bam_file} \
        -O raw.vcf.gz \

        
    


    
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

    a_write_log(log_text)
 
    index_ch = b_INDEXING(params.reference)
  
    def detected_contig = "chrM"

    c_MAPPING_2_SAM(params.reference, index_ch, read_pairs_ch)
    mapping_ch = c_MAPPING_2_SAM.out.mapping_test
    
    cleaned_ch = d_REMOVE_SOFT_CLIPPED_BASES(mapping_ch, params.python_script)
    // cleaned_ch.waitFor()

    fastq_ch = e_BACK_2_FASTQ(cleaned_ch)
    // fastq_ch.waitFor()

    merging_ch = f_MERGING(fastq_ch)
    // merging_ch.waitFor()

    trimming_ch = g_TRIMMING(merging_ch, params.left_primers, params.right_primers_rc)
    // trimming_ch.waitFor()

    mapping_final_ch = ha_MAPPING_2_BAM_trimmed(params.reference, index_ch, trimming_ch, params.amplicon_middle_positions)
    // mapping_final_ch.waitFor()

    mapping2_final_ch = ha2_MAPPING_2_BAM_merged(params.reference, index_ch, merging_ch, params.amplicon_middle_positions)
    // mapping_final_ch.waitFor()


    rtn_ch = hb_NUMTs_trimmed(mapping_final_ch, params.amplicon_middle_positions, params.humans, params.humans_amb, params.humans_ann, params.humans_bwt, params.humans_pac, params.humans_sa, params.numts, params.numts_amb, params.numts_ann, params.numts_bwt, params.numts_pac, params.numts_sa)
    // rtn_ch.waitFor()

    rtn2_ch = hb2_NUMTs_merged(mapping2_final_ch, params.amplicon_middle_positions, params.humans, params.humans_amb, params.humans_ann, params.humans_bwt, params.humans_pac, params.humans_sa, params.numts, params.numts_amb, params.numts_ann, params.numts_bwt, params.numts_pac, params.numts_sa)
    // rtn_ch.waitFor()


    i_ch = i_CALCULATE_STATISTICS(rtn_ch, params.python_script3)
    // i_ch.waitFor()

    j_ch = j_INDEX_CREATION(params.reference, detected_contig )
    // j_ch.waitFor()

    k_ch = k_MUTECT2(rtn_ch, j_INDEX_CREATION.out.ref_ch, j_INDEX_CREATION.out.fasta_index_ch, detected_contig, "mutect2_fusion")
    // k_ch.waitFor()

    vcf_ch = k_MUTECT2.out.mutect2_ch 

    // l_FINAL_VARIANTS(vcf_ch, params.reference, params.python_script2)
    
    fdstools_ch = FDSTOOLS(rtn2_ch)

    final_inputs = vcf_ch.join(fdstools_ch, by: 0)
    //     .map { sid, vcf_parts, fdstools_parts -> tuple(sid, *vcf_parts, *fdstools_parts) }

    // final_inputs = vcf_ch
    // .join(fdstools_ch, by: 0)
    // .map { sid, vcf_parts, fdstools_parts ->
    //     tuple(sid, *vcf_parts, *fdstools_parts)  // <-- this flattens the list properly
    // }

    // final_inputs.view()
    // final_inputs.into { final_variants_inputs }

    l_ch = l_FINAL_VARIANTS(final_inputs, params.reference, params.python_script2, params.python_script4, params.python_script5)
    variant_callers_ch = l_FINAL_VARIANTS.out.fdstools_mutect2_variants_ch

    n_ch = n_MERGE_FDSTOOLS_MUTECT2(variant_callers_ch, params.python_script6)

}