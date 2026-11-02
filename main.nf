nextflow.enable.dsl=2

params.reads = "$baseDir/raw_data/*_{R1,R2}_001.fastq.gz"
params.reference = "$baseDir/resources/rCRS/rCRS_NimaGen.fasta"

// flash
params.min_overlap = 10 // default in FLASH is 10
params.max_overlap = 226 // the read length: FLASH counts mismatches over at most this many overlap bases (default 65), and no overlap is longer than a read
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
params.minimum_length = 60
params.maximum_length = 600 // longer than any merged 2x226 read, so products spanning two amplicons are kept
params.primer_min_overlap = 15 // bases of a primer that must be in the read when it is incomplete at the read end

// bwa mem in p07 (primer-trimmed merged reads, judged by rtn in p09)
params.clipping_penalty = "100,100" // bwa mem -L; bwa default is 5,5

params.humans_index_dir = "$baseDir/resources/rtn_files/humans"
params.humans_index_base = "humans_NimaGen_2026-10-19.fa" // built by resources/rtn_files/humans/build/ from public sources

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
 
params.baseQ = 30 // minimum base quality: FDSTOOLS (TSSV) distrusts lower bases

params.python_script_remove_scb = "$baseDir/resources/scripts/remove_soft_clipped_bases_improved.py"
params.python_script_generate_read_depth_plot = "$baseDir/resources/scripts/generate_read_depth_plot_improved.py"
params.python_script_process_fdstools_sast = "$baseDir/resources/scripts/process_fdstools_output_improved_better.py"
params.python_script_merge_fdstools_mutect2 = "$baseDir/resources/scripts/merge_fdstools_mutect2_improved.py"
params.python_script_rtn_qualities = "$baseDir/resources/scripts/rtn_qualities.py"

// merging
params.depth = 20
params.min_vf_FDS = 2
params.lh_thresh = 10 // symmetric: floor=10%, ceiling=1-10%=90% — below floor not reported, floor-ceiling reported as LHP (lowercase), above ceiling reported as major (uppercase)
params.frame = "separate" // separate | shared: how 16180-16193 and 300-315 are written (resources/scripts/frames.py)
params.variant_note = true // false: the report has no variant_note column

    // rm -r "$baseDir/work"
    // rm -r "$baseDir/results"
    // rm .nextflow.*

// Each step is a module in modules/ (one process per file). Included after the params above,
// so that every module sees them.
include { p00_pipeline_parameters } from './modules/p00_pipeline_parameters'
include { p01_index_reference_fasta } from './modules/p01_index_reference_fasta'
include { p01b_prepare_humans_index } from './modules/p01b_prepare_humans_index'
include { p02_map_raw_fastq_p01 } from './modules/p02_map_raw_fastq_p01'
include { p03_filter_softclipped_fastq_p01_p02 } from './modules/p03_filter_softclipped_fastq_p01_p02'
include { p04_merge_fastq_p03 } from './modules/p04_merge_fastq_p03'
include { p05_trim_merged_fastq_p04 } from './modules/p05_trim_merged_fastq_p04'
include { p06_map_merged_bam_p01_p04 } from './modules/p06_map_merged_bam_p01_p04'
include { p07_map_merged_trimmed_bam_p01_p05 } from './modules/p07_map_merged_trimmed_bam_p01_p05'
include { p08_filter_numts_merged_fastq_p06 } from './modules/p08_filter_numts_merged_fastq_p06'
include { p09_filter_numts_trimmed_merged_bam_p07 } from './modules/p09_filter_numts_trimmed_merged_bam_p07'
include { p10_quality_control_p09 } from './modules/p10_quality_control_p09'
include { p11_variant_calling_fdstools_sast_p08 } from './modules/p11_variant_calling_fdstools_sast_p08'
include { p13_merge_variants_p10_p11 } from './modules/p13_merge_variants_p10_p11'

// A path inside the pipeline folder, written relative to it, so the parameter log stays short
def relativePath(path) {
    return path.toString().replace("${projectDir}/", "")
}

workflow {
    log_text = """\
         m t D N A - N i m a G e n  P I P E L I N E    
         ==========================================
         pipeline version                 : ${workflow.manifest.version}
         --reference                      : ${relativePath(params.reference)} # mtDNA reference genome
         --reads                          : ${relativePath(params.reads)}
         
         MERGING (with FLASH)
         --min_overlap                    : $params.min_overlap # The minimum required overlap length between two reads to provide a confident overlap (default: 10bp) 
         --max_overlap                    : $params.max_overlap # Mismatches counted over at most this many overlap bases; set to the read length
         --max_mismatch_density           : $params.max_mismatch_density # Maximum allowed ratio between the number of mismatched base pairs and the overlap length 

         TRIMMING (with CUTADAPT) 
         --adapter                        : $params.adapter # Anchor before the poly-C tract (278-302); R1 is cut after it (p02)
         --left_primers                   : ${relativePath(params.left_primers)} # Left primers (soft clips p03, trimming p05)
         --right_primers_rc               : ${relativePath(params.right_primers_rc)} # Right primers, reverse complement (p03, p05)
         --primer_min_overlap             : $params.primer_min_overlap # Bases of a primer that must be in the read when it is incomplete at the read end.
         --minimum_length                 : $params.minimum_length # Reads shorter than this are discarded
         --maximum_length                 : $params.maximum_length # Reads longer than this are discarded
         --discard-untrimmed              : enabled (hard coded) # Discard reads that do not contain an adapter/primer.

         MAPPING (with BWA MEM, p07)
         --clipping_penalty               : $params.clipping_penalty # Clipping penalty (-L) for the primer-trimmed merged reads rtn judges (bwa default: 5,5)
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
         --fdstools_library               : ${relativePath(params.fdstools_library)} # Position and flanks of each amplicon (p11; amplicon ranges p13)
         (tssv)
         --minimum                        : $params.minimum # Reads a sequence needs
         --num_threads                    : $params.num_threads # Threads of FDSTOOLS tssv, at most the cores of the computer
         --baseQ                          : $params.baseQ # Minimum base quality: FDSTOOLS distrusts lower bases

         POSTPROCESSING
         --depth                          : $params.depth # Read depth below which an amplicon is reported as LOW; its calls are kept
         --min_vf_FDS                     : $params.min_vf_FDS # Minor variant frequency threshold FDSTOOLS
         --lh_thresh                      : $params.lh_thresh # Length variants: under it not reported, up to 100 minus it lowercase, above uppercase
         --frame                          : $params.frame # separate (run lengths) or shared (laid from the left) for 16180-16193 and 300-315
         --variant_note                   : $params.variant_note # false: the report has no variant_note column

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


    // ────────────── PROCESS & MERGE VARIANTS ────────────────
    p13_merge_variants_p10_p11(p10_fdstools_ch, params.reference, params.python_script_process_fdstools_sast, params.python_script_merge_fdstools_mutect2,
                               files("$baseDir/resources/scripts/{notation,frames,iupac}.py"))

}
