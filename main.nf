nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/FASTQ/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/data/rCRS.fasta"
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 140 // default in FLASH is 65
params.max_mismatch_density = 0.25 //default in FLASH is 0.25
// params.multiqc = "$baseDir/multiqc"
params.publish_dir_mode = "symlink"
params.outdir = "results"

params.left_primers = "$baseDir/primers/left_primers.fasta"
params.right_primers_rc = "$baseDir/primers/right_primers_rc.fasta"


params.quality_cutoff = 25
params.minimum_length = 60
params.maximum_length = 300

params.detection_limit = 0.02
params.mapQ = 20
params.baseQ = 20
params.python_script = "$baseDir/scripts/remove_soft_clipped_bases.py"


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

         MUTSERVE 
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
process write_log{
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

process INDEX {
    tag "bwa index on $reference"   
    publishDir "$params.outdir/index", mode: 'copy'

    input:
    path reference

    output:
    path("${reference}.*")

    script:
    """
    bwa index $reference 
    """
}

process MAPPING_2_SAM {
    tag "bwa mem on $sample_id"
    publishDir "$params.outdir/mapped_sam", mode: 'copy'

    input:
    path reference
    path index_files
    tuple val(sample_id), path(reads)
    // path merged_file
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.sam"), path("${sample_id}_R2.sam")

    script:
    """
    bwa mem $reference ${reads[0]} > ${sample_id}_R1.sam 
    bwa mem $reference ${reads[1]} > ${sample_id}_R2.sam 
    """
}

process REMOVE_SOFT_CLIPPED_BASES {
    tag "remove_scb on $sample_id"
    publishDir "$params.outdir/cleaned", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_r1), path(sam_r2)
    // path merged_file
    
    output:
    tuple val(sample_id), path("${sample_id}_R1_cleaned.bam"), path("${sample_id}_R2_cleaned.bam")


    script:
    """
    cat ${sam_r1} | python $params.python_script > ${sample_id}_R1_cleaned.sam 
    cat ${sam_r2} | python $params.python_script > ${sample_id}_R2_cleaned.sam
    samtools-1.21 view -Sb ${sample_id}_R1_cleaned.sam >  ${sample_id}_R1_cleaned.bam
    samtools-1.21 view -Sb ${sample_id}_R2_cleaned.sam >  ${sample_id}_R2_cleaned.bam
    """
}

process BACK_2_FASTQ {
    tag "convert_2_fastq on $sample_id"
    publishDir "$params.outdir/fastq", mode: 'copy'

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




process MERGING {
    tag "flash on $sample_id"
    publishDir "$params.outdir/merged", mode: 'copy'

    input:
    tuple val(sample_id), path(cleaned_fastq_r1), path(cleaned_fastq_r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_cleaned_merged.extendedFrags.fastq")

    script:
    """
    flash ${cleaned_fastq_r1} ${cleaned_fastq_r2} -m $params.min_overlap -M $params.max_overlap -x $params.max_mismatch_density -O -o ${sample_id}_cleaned_merged
    """
}

process TRIMMING {
    tag "cutadapt on $sample_id"
    publishDir "$params.outdir/cutadapt", mode: 'copy'

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

process MAPPING_2_BAM {
    tag "bwa mem on $sample_id"
    publishDir "$params.outdir/mapped_final", mode: 'copy'

    input:
    path reference
    path index_files
    tuple val(sample_id), path(trimmed_fastq)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")

    script:
    """
    bwa mem $reference ${trimmed_fastq} | samtools-1.21 view -Sb - > ${sample_id}_tmp.bam    
    samtools-1.21 view -h ${sample_id}_tmp.bam | awk '\$1 ~ /^@/ || \$6 !~ /S/' | samtools-1.21 view -b -o ${sample_id}.bam
    samtools-1.21 sort -o ${sample_id}_sorted.bam ${sample_id}.bam
    samtools-1.21 index ${sample_id}_sorted.bam
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
        --strand-bias 1.6 \
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


// process MUTECT2 {

//     input:
//     path bam_file
//     path reference
//     // path fasta_index_files
//     val method

//     output:
//     tuple path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), val(method), emit: mutect2_ch

//     script:
//     def avail_mem = 1024
//     if (task.memory) {
//         avail_mem = (task.memory.mega*0.8).intValue()
//     }    

//     """
//     samtools-1.21 index ${bam_file}
//     samtools-1.21 faidx ${reference}
//     samtools-1.21 dict ${reference} -o rCRS.dict
//     gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \
//         Mutect2 \
//         -R ${reference} \
//         --min-base-quality-score ${params.baseQ} \
//         -callable-depth 6 \
//         --native-pair-hmm-threads 6 \
//         --max-reads-per-alignment-start 0 \
//         --tmp-dir . \
//         -I ${bam_file} \
//         -O raw.vcf.gz
    
//     gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \
//         FilterMutectCalls \
//         -R ${reference} \
//         --min-reads-per-strand 2 \
//         -V raw.vcf.gz \
//         --tmp-dir . \
//         -O ${bam_file.baseName}.vcf.gz

//     bcftools norm \
//         -m-any \
//         -f ${reference} \
//         -o ${bam_file.baseName}.norm.vcf.gz -Oz \
//         ${bam_file.baseName}.vcf.gz 

//     bcftools view \
//     -i 'FORMAT/AF>=${params.detection_limit}' \
//     -o ${bam_file.baseName}.vcf.gz -Oz \
//     ${bam_file.baseName}.norm.vcf.gz 
    
//     tabix -f ${bam_file.baseName}.vcf.gz

//     rm ${bam_file.baseName}.norm.vcf.gz 
//     rm raw.vcf.gz
//     """
// }


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    write_log(log_text)
    index_ch = INDEX(params.reference)
    // merging_ch = MERGING(read_pairs_ch)
    mapping_ch = MAPPING_2_SAM(params.reference, index_ch, read_pairs_ch)
    cleaned_ch = REMOVE_SOFT_CLIPPED_BASES(mapping_ch)
    fastq_ch = BACK_2_FASTQ(cleaned_ch)
    merging_ch = MERGING(fastq_ch)
    trimming_ch = TRIMMING(merging_ch)
    mapping_final_ch = MAPPING_2_BAM(params.reference, index_ch, trimming_ch)
    // INDEX_CREATION(params.reference)
    MUTSERVE(params.reference, index_ch, "mutserve_fusion", mapping_final_ch)

    // MUTECT2(mapping_ch, params.reference, "mutect2_fusion")
        
    // vcf_ch = MUTSERVE.out.mutserve_ch.concat(MUTECT2.out.mutect2_ch)
    // file_count =  MUTSERVE.out.mutserve_ch.count()

    
}