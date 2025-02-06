nextflow.enable.dsl=2

/* 
 * pipeline input parameters 
 */
params.reads = "$baseDir/data/FASTQ/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/data/rCRS.fasta"
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 65 // default in FLASH is 65
params.max_mismatch_density = 0.25 //default in FLASH is 0.25
// params.multiqc = "$baseDir/multiqc"
params.publish_dir_mode = "symlink"
params.outdir = "results"

params.detection_limit = 0.02
params.mapQ = 20
params.baseQ = 20

    // rm -r "$baseDir/work"
    // rm -r "$baseDir/results"
    // rm .nextflow.*

log_text = """\
         m t D N A - N i m a G e n  P I P E L I N E    
         ==========================================
         mtDNA reference genome : ${params.reference}
         reads                  : ${params.reads}
         
         MERGING (with FLASH)
         merging_min-overlap    : $params.min_overlap
         merging_max-overlap    : $params.max_overlap
         max-mismatch-density   : $params.max_mismatch_density
         
         MUTSERVE 
         detection_limit        : $params.detection_limit
         mapQ                   : $params.mapQ
         baseQ                  : $params.baseQ

         publish_dir_mode       : $params.publish_dir_mode
         outdir                 : ${params.outdir}
         """

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
    publishDir "$params.outdir/index", mode: params.publish_dir_mode

    input:
    path reference

    output:
    path("${reference}.*")

    script:
    """
    bwa index $reference 
    """
}

process MERGING {
    tag "flash on $sample_id"
    publishDir "$params.outdir/merged", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(reads) // Paired-end FASTQ files
    
    output:
    path "${sample_id}_merged.extendedFrags.fastq"

    script:
    """
    flash ${reads[0]} ${reads[1]} -M $params.max_overlap -x $params.max_mismatch_density -o ${sample_id}_merged
    """
}

process MAPPING {
    tag "bwa mem on ${merged_file[0].baseName}"
    publishDir "$params.outdir/mapped", mode: params.publish_dir_mode

    input:
    path reference
    path index_files
    path merged_file
    
    output:
    path "${merged_file[0].baseName}_sorted.bam"
    path "${merged_file[0].baseName}_sorted.bai"

    script:
    """
    bwa mem ${index_files[0].baseName} ${merged_file} | samtools-1.21 view -Sb - > ${merged_file[0].baseName}.bam
    samtools-1.21 sort -o ${merged_file[0].baseName}_sorted.bam ${merged_file[0].baseName}.bam
    samtools-1.21 index ${merged_file[0].baseName}_sorted.bam
    """
    // bwa mem ${index_files[0].baseName} ${merged_file} | samtools sort -o ${merged_file[0].baseName}.bam
    // samtools-1.21 view -h ${merged_file[0].baseName}_sorted.bam | '\\\$1 ~ /^@/ || \\\$6 !~ /S/' | samtools-1.21 view -b -o ${merged_file[0].baseName}__sorted_filtered.bam  
    // samtools index ${merged_file[0].baseName}.bam
    // samtools-1.21 view -h -F 4 ${merged_file[0].baseName}_sorted.bam | grep -v 'S' | samtools-1.21 view -b > ${merged_file[0].baseName}_sorted_filtered.bam  
    // samtools-1.21 index ${merged_file[0].baseName}_sorted_filtered.bam
}

// process MUTSERVE {

//     input:
//     path bam_file
//     path reference
//     val method

//     output:
//     tuple path("${bam_file.baseName}.vcf.gz"), path("${bam_file.baseName}.vcf.gz.tbi"), val(method), emit: mutserve_ch
    
//     script:
//     def avail_mem = 1024
//     if (task.memory) {
//         avail_mem = (task.memory.mega*0.8).intValue()
//     }    

//     """
//     #todo: check used mutserve strand-bias with default parameter
//     samtools-1.21 index ${bam_file} 
//     java -Xmx${avail_mem}M -jar /opt/mutserve/mutserve.jar \
//         call \
//         --level ${params.detection_limit} \
//         --reference ${reference} \
//         --mapQ ${params.mapQ} \
//         --baseQ ${params.baseQ} \
//         --output ${bam_file.baseName}.vcf.gz \
//         --no-ansi \
//         --strand-bias 1.6 \
//         --write-raw \
//         ${bam_file} 

//     bcftools norm \
//         -m-any \
//         -f ${reference} \
//         -o ${bam_file.baseName}.norm.vcf.gz -Oz \
//         ${bam_file.baseName}.vcf.gz 
    
//     mv ${bam_file.baseName}.norm.vcf.gz ${bam_file.baseName}.vcf.gz
//     tabix -f ${bam_file.baseName}.vcf.gz
//     """
// }

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
    merging_ch = MERGING(read_pairs_ch)
    mapping_ch = MAPPING(params.reference, index_ch, merging_ch)
    
    // INDEX_CREATION(params.reference)
    // MUTSERVE(mapping_ch, params.reference, "mutserve_fusion")

    // MUTECT2(mapping_ch, params.reference, "mutect2_fusion")
        
    // vcf_ch = MUTSERVE.out.mutserve_ch.concat(MUTECT2.out.mutect2_ch)
    // file_count =  MUTSERVE.out.mutserve_ch.count()

    
}