nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/FASTQ/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/data/rCRS2.fasta"
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 140 // default in FLASH is 65
params.max_mismatch_density = 0.25 //default in FLASH is 0.25
// params.multiqc = "$baseDir/multiqc"
params.publish_dir_mode = "symlink"
params.outdir = "results"

params.adapter = 'ATCATAACAAAAAATTTCCACCAAA'

params.left_primers = "$baseDir/primers/left_primers.fasta"
params.right_primers_rc = "$baseDir/primers/right_primers_rc.fasta"
params.amplicon_middle_positions = "$baseDir/primers/amplicons_bed.txt"

// cutadapt
params.quality_cutoff = 25
params.minimum_length = 60
params.maximum_length = 300

params.humans = "/Users/peter/Documents/Postdoc/Projects/for_Mario/mtDNA_NimaGen/RtN/humans.fa"
params.numts = "/Users/peter/Documents/Postdoc/Projects/for_Mario/mtDNA_NimaGen/RtN/Calabrese_Dayama_Smart_Numts.fa"


// mutet2
params.detection_limit = 0.075
params.mapQ = 30
params.baseQ = 32
params.alignQ = 30
params.mode = 'fusion'

params.python_script = "$baseDir/scripts/remove_soft_clipped_bases.py"
params.python_script2 = "$baseDir/scripts/python_empop.py"
params.python_script3 = "$baseDir/scripts/python_coverage.py"

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

process b_INDEX {
    tag "b: bwa index on $reference"   
    publishDir "$params.outdir/b_index", mode: 'copy'

    input:
    path reference

    output:
    path("${reference}.*")

    script:
    """
    bwa index $reference 
    """
}

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
    // mv ${reads[0]} tmp.fastq.gz
    // cutadapt -a $params.adapter -o ${reads[0]}  tmp.fastq.gz 

process d_REMOVE_SOFT_CLIPPED_BASES {
    container 'peterresutik/nimagen-pipeline:latest'
    tag "d: remove_scb on $sample_id"
    publishDir "$params.outdir/d_cleaned", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_r1), path(sam_r2)
    // path merged_file
    
    output:
    tuple val(sample_id), path("${sample_id}_R1_cleaned.bam"), path("${sample_id}_R2_cleaned.bam")
    

    script:
    """
    cat ${sam_r1} | python /workspace/scripts/remove_soft_clipped_bases.py > ${sample_id}_R1_cleaned.sam 
    cat ${sam_r2} | python /workspace/scripts/remove_soft_clipped_bases.py > ${sample_id}_R2_cleaned.sam
    samtools view -Sb ${sample_id}_R1_cleaned.sam >  ${sample_id}_R1_cleaned.bam
    samtools view -Sb ${sample_id}_R2_cleaned.sam >  ${sample_id}_R2_cleaned.bam
    """
}

process e_BACK_2_FASTQ {
    tag "e: convert_2_fastq on $sample_id"
    publishDir "$params.outdir/e_fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(cleaned_bam_r1), path(cleaned_bam_r2)
    // path merged_file
    
    output:
     tuple val(sample_id), path("${sample_id}_R1_cleaned.fastq"), path("${sample_id}_R2_cleaned.fastq")

    script:
    """
    samtools-1.21 fastq ${cleaned_bam_r1} > ${sample_id}_R1_cleaned.fastq    
    samtools-1.21 fastq ${cleaned_bam_r2} > ${sample_id}_R2_cleaned.fastq   
    """
}




process f_MERGING {
    tag "f: flash on $sample_id"
    publishDir "$params.outdir/f_merged", mode: 'copy'

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
    publishDir "$params.outdir/g_cutadapt", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_cleaned_merged_trimmed_left_right.fastq")

    script:
    """
    cutadapt -g file:$params.left_primers -q $params.quality_cutoff -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_cleaned_merged_trimmed_left.fastq $merged_fastq  
    cutadapt -a file:$params.right_primers_rc -q $params.quality_cutoff -m $params.minimum_length -M $params.maximum_length --discard-untrimmed -o ${sample_id}_cleaned_merged_trimmed_left_right.fastq ${sample_id}_cleaned_merged_trimmed_left.fastq    
    """
}
// --discard-untrimmed

process ha_MAPPING_2_BAM {
    tag "ha: bwa mem on $sample_id"
    publishDir "$params.outdir/h_mapped_final", mode: 'copy'

    input:
    path reference
    path index_files
    tuple val(sample_id), path(trimmed_fastq)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}_coverage.txt")

    script:
    """
    bwa mem $reference ${trimmed_fastq} | samtools-1.21 view -Sb - > ${sample_id}_tmp.bam    
    samtools-1.21 view -h ${sample_id}_tmp.bam | awk '\$1 ~ /^@/ || \$6 !~ /S/' | samtools-1.21 view -b -o ${sample_id}.bam
    samtools-1.21 sort -o ${sample_id}_tmp.bam ${sample_id}.bam
    
    samtools-1.21 addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${sample_id}_tmp.bam  -o ${sample_id}_tmp2.bam
    mv ${sample_id}_tmp2.bam ${sample_id}.bam
    samtools-1.21 index ${sample_id}.bam

    samtools-1.21 depth -a -b $params.amplicon_middle_positions ${sample_id}.bam > ${sample_id}_coverage.txt
    
    """
}

process hb_NUMTs {
    tag "hb: rnt on $sample_id"
    publishDir "$params.outdir/hb_numts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt)

    output:
    tuple val(sample_id), path("${sample_id}.rtn.bam"), path("${sample_id}.rtn.bam.bai"), path("${sample_id}_coverage_numts.txt")

    script:
    """
    /Users/peter/Documents/Postdoc/Projects/for_Mario/mtDNA_NimaGen/RtN/rtn -p -h $params.humans -n $params.numts -b $bam_file



    samtools-1.21 view -h -q 30 ${sample_id}.rtn.bam > ${sample_id}.rtn_tmp.bam

    samtools-1.21 depth -a -b $params.amplicon_middle_positions ${sample_id}.rtn_tmp.bam > ${sample_id}_coverage_numts.txt
    
    """
}

// process 2nd_MAPPING {
//     tag "bwa mem on ${merged_file[0].baseName}"
//     publishDir "$params.outdir/mapped", mode: params.publish_dir_mode

//     input:
//     path reference
//     path index_files
//     path merged_file
    
//     output:
//     path "${merged_file[0].baseName}_sorted.bam"
//     path "${merged_file[0].baseName}_sorted.bam.bai"

//     script:
//     """
//     bwa mem ${index_files[0].baseName} ${merged_file} | samtools-1.21 view -Sb - > ${merged_file.baseName}.bam
//     samtools-1.21 sort -o ${merged_file[0].baseName}_sorted.bam ${merged_file[0].baseName}.bam
//     samtools-1.21 index ${merged_file[0].baseName}_sorted.bam
//     """
//     // | samtools view -Sb - > ${merged_file.baseName}.bam
//     // samtools sort -o ${merged_file.baseName}_sorted.bam ${merged_file.baseName}.bam
//     // samtools index ${merged_file.baseName}_sorted.bam
//     // bwa mem ${index_files[0].baseName} ${merged_file} | samtools-1.21 view -Sb - > ${merged_file.baseName}.bam
//     // bwa mem ${index_files[0].baseName} ${merged_file} | samtools-1.21 view -Sb - > ${merged_file[0].baseName}.bam
//     // samtools-1.21 sort -o ${merged_file[0].baseName}_sorted.bam ${merged_file[0].baseName}.bam
//     // samtools-1.21 index ${merged_file[0].baseName}_sorted.bam
    
//     // bwa mem ${index_files[0].baseName} ${merged_file} | samtools sort -o ${merged_file[0].baseName}.bam
//     // samtools-1.21 view -h ${merged_file[0].baseName}_sorted.bam | '\\\$1 ~ /^@/ || \\\$6 !~ /S/' | samtools-1.21 view -b -o ${merged_file[0].baseName}__sorted_filtered.bam  
//     // samtools index ${merged_file[0].baseName}.bam
//     // samtools-1.21 view -h -F 4 ${merged_file[0].baseName}_sorted.bam | grep -v 'S' | samtools-1.21 view -b > ${merged_file[0].baseName}_sorted_filtered.bam  
//     // samtools-1.21 index ${merged_file[0].baseName}_sorted_filtered.bam
// }

process MUTSERVE {
    tag "mutserve on $sample_id"
    publishDir "$params.outdir/mutserve", mode: 'copy'

    input:
    path reference
    path index_files
    val method
    tuple val(sample_id), path(sorted_bam), path(sorted_index)

    output:
    tuple path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), val(method), emit: mutserve_ch
    
    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    

    """
    #todo: check used mutserve strand-bias with default parameter 
    java -Xmx${avail_mem}M -jar /opt/mutserve/mutserve.jar \
        call \
        --level ${params.detection_limit} \
        --reference ${reference} \
        --mapQ ${params.mapQ} \
        --baseQ ${params.baseQ} \
        --output ${sample_id}.vcf.gz \
        --no-ansi \
        --write-raw \
        ${sorted_bam} 

    bcftools norm \
        -m-any \
        -f ${reference} \
        -o ${sample_id}.norm.vcf.gz -Oz \
        ${sample_id}.vcf.gz 
    
    mv ${sample_id}.norm.vcf.gz ${sample_id}.vcf.gz
    tabix -f ${sample_id}.vcf.gz
    """
}

// process INDEX_CREATION {
// 	input:
// 	path reference
    
// 	output:
// 	path "rCRS*.{dict,fai}", emit: fasta_index_ch
//     // path "ref*.{dict,fai}", emit: fasta_index_ch
// 	// path "ref.fasta", emit: ref_ch
	
//     """
//     samtools-1.21 faidx rCRS.fasta
//     samtools-1.21 dict rCRS.fasta -o rCRS.dict
// 	"""
// }

process i_CALCULATE_STATISTICS {
    tag "i: calculate_statistics on $sample_id"
    publishDir "$params.outdir/i_calculate_statistics", mode: 'copy'

    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt)
    tuple val(sample_id2), path(bam_file2), path(bam_index2), path(coverage_numts_txt)

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
    // samtools-1.21 addreplacerg -r '@RG\tID:${sample_id}\tSM:${sample_id}' ${bam_file}  -o ${sample_id}_tmp.bam
    // mv ${sample_id}_tmp.bam ${bam_file}

    """
    ## Create Mapping File
    ## mkdir "${sample_id}"
    ## echo -e "Sample\tFilename" > /${sample_id}/$mapping_name
    echo -e "Sample\tFilename" > $mapping_name

    echo "\$(samtools-1.21 samples ${bam_file})" >> $mapping_name

    ## Calculate summary statistics
    samtools-1.21 coverage ${bam_file} > samtools_coverage_${sample_id}.txt
    csvtk grep -t -f3 -p 16623 -C '\$' samtools_coverage_${sample_id}.txt -T -o mtdna.txt --num-cpus ${task.cpus} 
        
    contig=\$(csvtk cut -t -f 1 mtdna.txt --num-cpus ${task.cpus})
    numreads=\$(csvtk cut -t -f 4 mtdna.txt --num-cpus ${task.cpus})
    covered_bases=\$(csvtk cut -t -f 5 mtdna.txt --num-cpus ${task.cpus})
    covered_bases_percentage=\$(csvtk cut -t -f 6 mtdna.txt --num-cpus ${task.cpus})
    mean_depth=\$(csvtk cut -t -f 7 mtdna.txt --num-cpus ${task.cpus})
    mean_base_quality=\$(csvtk cut -t -f 8 mtdna.txt --num-cpus ${task.cpus})
    mean_map_quality=\$(csvtk cut -t -f 9 mtdna.txt --num-cpus ${task.cpus})
    readgroup=\$(samtools-1.21 view -H ${bam_file} | csvtk grep -I -H -r -p "^@RG" --num-cpus ${task.cpus} | sed 's/\t/,/g' | head -n 1)
    
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

    python $params.python_script3 $coverage_txt $coverage_numts_txt ${sample_id}_coverage_plot.png

    """
}
// echo "\$(${sample_id} ${bam_file})" >> $mapping_name

process INPUT_VALIDATION {
    // tag "input_validation on $sample_id"
    publishDir "$params.outdir/input_validation", mode: 'copy'

    input:
    path bams_ch
    path statistics
    path mapping
    path reference
    path index_files

    output:
    path("sample_statistics.txt"), emit: summarized_ch
    path("sample_mappings.txt"), emit: mapping_ch
    path("excluded_samples.txt"), emit: excluded_ch
    path("contig.txt"), emit: contig_ch
    path("*.bam"), includeInputs: true, emit: validated_files

    """
    csvtk concat \
        -t ${statistics} \
        -T -o sample_statistics.txt \
        --num-cpus ${task.cpus}
    
    csvtk concat \
        -t ${mapping} \
        -T -o sample_mappings.txt \
        --num-cpus ${task.cpus}
    
    java -jar /opt/mutserve/mutserve.jar stats \
        --input sample_statistics.txt \
        --mapping sample_mappings.txt \
        --detection-limit ${params.detection_limit}  \
        --reference ${reference}  \
        --baseQ ${params.baseQ}\
        --mapQ ${params.mapQ} \
        --alignQ ${params.alignQ} \
        --output-excluded-samples excluded_samples.txt \
        --output-contig contig.txt \
        --tool ${params.mode}

    # delete excluded_samples from BAM input channel directly
    awk -v q='"' '{print "rm " q \$1 q }' excluded_samples.txt | sh

   
    """
}
//  python -m json.tool cloudgene.report.json




// process QUALITY_CONTROL {
    
//     publishDir "${params.output_reports}/multiqc", mode: "copy", pattern: '*.html'
    
//     input:
//     path zip
    
//     output:
// 	path "*.html"

// 	"""
// 	multiqc . --filename index.html
// 	"""
// }

process j_INDEX_CREATION {
	publishDir "$params.outdir/j_index_creation", mode: 'copy'
    input:
	path reference
	val mtdna_tag

	output:
	path "ref*.{dict,fai}", emit: fasta_index_ch
	path "ref.fasta", emit: ref_ch

	"""
	sed -e "s/^>.*/>$mtdna_tag/" $reference > ref.fasta
    samtools-1.21 faidx ref.fasta
    	samtools-1.21 dict ref.fasta \
	    -o ref.dict
	"""
}


process k_MUTECT2 {
    tag "k: mutect2 on $sample_id"
    publishDir "$params.outdir/k_mutect2", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index), path(coverage_txt)
    path reference
    path fasta_index_files
    val detected_contig
    val method

    output:
    tuple  val(sample_id), path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), val(method), emit: mutect2_ch

    script:
    def avail_mem = 1024
    if (task.memory) {
        avail_mem = (task.memory.mega*0.8).intValue()
    }    
    // samtools-1.21 index ${bam_file}
    
    """

    /Users/peter/anaconda3/pkgs/gatk4-4.6.1.0-py310hdfd78af_0/share/gatk4-4.6.1.0-0/gatk  --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \
        Mutect2 \
        -R ${reference} \
        -L '${detected_contig}' \
        --min-base-quality-score ${params.baseQ} \
        -callable-depth 6 \
        --native-pair-hmm-threads ${task.cpus} \
        --max-reads-per-alignment-start 0 \
        --tmp-dir . \
        -I ${bam_file} \
        -O raw.vcf.gz
    
    /Users/peter/anaconda3/pkgs/gatk4-4.6.1.0-py310hdfd78af_0/share/gatk4-4.6.1.0-0/gatk  --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \
        FilterMutectCalls \
        -R ${reference} \
        --min-reads-per-strand 2 \
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
    rm raw.vcf.gz
    """
}

process l_ANNOTATE_VARIANTS {
    tag "l: annotate_variants on $sample_id"
    publishDir "$params.outdir/l_filter_variants", mode: 'copy'
    // publishDir "${params.output}", mode: 'copy'

    input:
    tuple  val(sample_id), path(vcf_file), path(vcf_file_idx), val(method)
    path reference

    output:
    path("${vcf_file.baseName}.${method}.filtered.empop.txt"), emit: combined_methods_ch

    script:
    def vcf_name = "${vcf_file}".replaceAll('.vcf.gz', '')

    """
    echo -e "test"
    echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tMeanBaseQuality\tCoverage\tGT" \
        > ${vcf_file.baseName}.${method}.txt

    bcftools query -u \
        -f '${vcf_name}.bam\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%BQ\t%DP\t%GT]\n' \
        ${vcf_file} >> ${vcf_file.baseName}.${method}.txt    
    

    ## annotating SNVS and INDELs for reporting
    awk 'BEGIN {OFS="\t"} {
        if (NR == 1) { print \$0, "Type"; next }
        if ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5)) { \$10="3" }
        else if (\$9 == "1") { \$10="1" }
        else if (\$9 == "0/1" || \$9 == "1/0" || \$9 == "0|1" || \$9 == "1|0") { \$10="2" }
        else { \$10="UNKNOWN" }
        print
    }' ${vcf_file.baseName}.${method}.txt > ${vcf_file.baseName}.${method}.filtered.txt

    python $params.python_script2 ${vcf_file.baseName}.${method}.filtered.txt ${vcf_file.baseName}.${method}.filtered.empop.txt $reference 
    
    """
}

    // if [[ ${method} == "mutserve_fusion" ]]
    // then
    //     awk -F'\t' 'NR == 1 || (length(\$4) == 1 && length(\$5) == 1)' \
    //         ${vcf_file.baseName}.${method}.txt > ${vcf_file.baseName}.${method}.filtered.tmp.txt

    // elif [[ ${method} == "mutect2_fusion" ]]
    // then
    //     awk -F'\t' 'NR == 1 || ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5))' \
    //         ${vcf_file.baseName}.${method}.txt > ${vcf_file.baseName}.${method}.filtered.tmp.txt
    // else 
    //     mv ${vcf_file.baseName}.${method}.txt ${vcf_file.baseName}.${method}.filtered.tmp.txt  
    // fi
    
    // ## annotating SNVS and INDELs for reporting
    // awk 'BEGIN {OFS="\t"} {
    //     if (NR == 1) { print \$0, "Type"; next }
    //     if ((length(\$4) > 1 || length(\$5) > 1) && length(\$4) != length(\$5)) { \$10="3" }
    //     else if (\$9 == "1") { \$10="1" }
    //     else if (\$9 == "0/1" || \$9 == "1/0" || \$9 == "0|1" || \$9 == "1|0") { \$10="2" }
    //     else { \$10="UNKNOWN" }
    //     print
    // }' ${vcf_file.baseName}.${method}.filtered.tmp.txt > ${vcf_file.baseName}.${method}.filtered.txt

    // rm ${vcf_file.baseName}.${method}.filtered.tmp.txt

process MERGING_VARIANTS {
    publishDir "$params.outdir/merged_variants", mode: 'copy'
    input:
    path variants_txt
    val mode

    output:
    path("variants.txt"), emit: txt_summarized_ch

    """
    csvtk concat \
        -t ${variants_txt} \
        -T -o variants.concat.txt \
        --num-cpus ${task.cpus}
    
    csvtk sort \
        -t variants.concat.txt \
        -k ID:N -k Pos:n -k Ref:N -k Type:nr  -k Variant:N \
        -T -o variants.sorted.txt \
        --num-cpus ${task.cpus}

    if [[ ${mode} == "fusion" ]]
    then
        java -cp "/opt/VariantMerger.jar:/opt/lib/*" VariantMerger \
            variants.sorted.txt \
            --output variants.txt
    else
        mv variants.sorted.txt variants.txt
    fi
    """
}

// process EMPOP {
//     publishDir "$params.outdir/empop_variants", mode: 'copy'
//     input:
//     path variants_txt
//     val mode

//     """
//     python empop_formatter.py variants.txt reference.fasta output.txt
//     """
// }


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    a_write_log(log_text)
    index_ch = b_INDEX(params.reference)
    def detected_contig = "chrM"
    // INPUT_VALIDATION.out.contig_ch.text.trim()


    // merging_ch = MERGING(read_pairs_ch)
    // mapping_ch = 
    c_MAPPING_2_SAM(params.reference, index_ch, read_pairs_ch)
    
    mapping_ch = c_MAPPING_2_SAM.out.mapping_test
    cleaned_ch = d_REMOVE_SOFT_CLIPPED_BASES(mapping_ch)
    // fastq_ch = e_BACK_2_FASTQ(cleaned_ch)
    // merging_ch = f_MERGING(fastq_ch)
    // trimming_ch = g_TRIMMING(merging_ch)
    // mapping_final_ch = ha_MAPPING_2_BAM(params.reference, index_ch, trimming_ch)

    // rnt_ch = hb_NUMTs(mapping_final_ch)

    // i_CALCULATE_STATISTICS(mapping_final_ch, rnt_ch)

    // // INPUT_VALIDATION(
    // //     CALCULATE_STATISTICS.out.fixed_file.collect(),
    // //     CALCULATE_STATISTICS.out.stats_ch.collect(),
    // //     CALCULATE_STATISTICS.out.mapping_ch.collect(),
    // //     params.reference, index_ch
    // // )



    // // QUALITY_CONTROL(
    // //     CALCULATE_STATISTICS.out.fastqc_ch.collect()
    // // )

    // // // haplogrep_ch = file("$projectDir/files/haplogroups.txt")
    // // // contamination_ch = file("$projectDir/files/haplocheck.txt")

    // // validated_files = INPUT_VALIDATION.out.validated_files.flatten()
     
    // // MUTSERVE(params.reference, index_ch, "mutserve_fusion", mapping_final_ch)

    // j_INDEX_CREATION(
    //     params.reference,
    //     detected_contig
    // )

    // k_MUTECT2(
    //     rnt_ch,
    //     j_INDEX_CREATION.out.ref_ch,
    //     j_INDEX_CREATION.out.fasta_index_ch,
    //     detected_contig,
    //     "mutect2_fusion"
    // )
    
    // // // vcf_ch = MUTSERVE.out.mutserve_ch
    // vcf_ch = k_MUTECT2.out.mutect2_ch 
    // // vcf_ch = MUTSERVE.out.mutserve_ch.concat(MUTECT2.out.mutect2_ch)
    // // file_count =  MUTSERVE.out.mutserve_ch.count()
    
    // l_ANNOTATE_VARIANTS (
    //     vcf_ch,
    //     params.reference
    // )

    // // MERGING_VARIANTS(
    // //     FILTER_VARIANTS.out.combined_methods_ch.collect(),
    // //     params.mode
    // // )
    
    // // variants_txt_ch = MERGING_VARIANTS.out.txt_summarized_ch    


    // // MUTECT2(mapping_ch, params.reference, "mutect2_fusion")
        
    // // vcf_ch = MUTSERVE.out.mutserve_ch.concat(MUTECT2.out.mutect2_ch)
    // // file_count =  MUTSERVE.out.mutserve_ch.count()

    
}