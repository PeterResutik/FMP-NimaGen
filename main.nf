nextflow.enable.dsl=2

// Each step is a module in modules/ (one process per file, numbered in the order the steps run); the parameters are in nextflow.config.
include { p00_parameters } from './modules/p00_parameters'
include { p01_index_reference } from './modules/p01_index_reference'
include { p02_index_humans } from './modules/p02_index_humans'
include { p03_map_read_pairs } from './modules/p03_map_read_pairs'
include { p04_remove_soft_clips } from './modules/p04_remove_soft_clips'
include { p05_merge_read_pairs } from './modules/p05_merge_read_pairs'
include { p06_trim_primers } from './modules/p06_trim_primers'
include { p07_map_merged_reads } from './modules/p07_map_merged_reads'
include { p08_map_trimmed_reads } from './modules/p08_map_trimmed_reads'
include { p09_filter_numts } from './modules/p09_filter_numts'
include { p10_reads_for_fdstools } from './modules/p10_reads_for_fdstools'
include { p11_read_depth } from './modules/p11_read_depth'
include { p12_fdstools } from './modules/p12_fdstools'
include { p13_report } from './modules/p13_report'

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
         --adapter                        : $params.adapter # Anchor before the poly-C tract (278-302); R1 is cut after it (p03)
         --left_primers                   : ${relativePath(params.left_primers)} # Left primers (soft clips p04, trimming p06)
         --right_primers_rc               : ${relativePath(params.right_primers_rc)} # Right primers, reverse complement (p04, p06)
         --primer_min_overlap             : $params.primer_min_overlap # Bases of a primer that must be in the read when it is incomplete at the read end.
         --minimum_length                 : $params.minimum_length # Reads shorter than this are discarded
         --maximum_length                 : $params.maximum_length # Reads longer than this are discarded
         --discard-untrimmed              : enabled (hard coded) # Discard reads that do not contain an adapter/primer.

         MAPPING (with BWA MEM, p08)
         --clipping_penalty               : $params.clipping_penalty # Clipping penalty (-L) for the primer-trimmed merged reads rtn judges (bwa default: 5,5)
         --amplicon_middle_positions      : ${relativePath(params.amplicon_middle_positions)} # Middle position of each amplicon, for its read depth

         NUMTs REMOVAL (with RTN)
         --humans_index_dir               : ${relativePath(params.humans_index_dir)} # Folder of the human mitogenomes rtn compares reads with
         --humans_index_base              : $params.humans_index_base # Human mitogenomes, unpacked and indexed on the first run
         --numts_index_dir                : ${relativePath(params.numts_index_dir)} # Folder of the NUMT sequences rtn compares reads with
         --numts_index_base               : $params.numts_index_base # NUMT sequences (from rtn)
         --rtn_positions                  : ${relativePath(params.rtn_positions)} # positions rtn sees at quality 2 (none = off)
         --rtn_likelihood                 : $params.rtn_likelihood # minimum likelihood of a read under the human mitogenomes (rtn -L; rtn default 1e-5)
         --mapQ                           : $params.mapQ # Used to filter out reads assigned as NUMTs by RTN 
         --skip_numt_filter               : $params.skip_numt_filter # true skips rtn and the end-to-end check (simulated reads contain no NUMTs); the MAPQ filter still applies
         --rtn_memory                     : $params.rtn_memory # memory per rtn task, which holds the humans index (about 2 GB)

         VARIANT CALLING (with FDSTOOLS)
         --fdstools_library               : ${relativePath(params.fdstools_library)} # Position and flanks of each amplicon (p12; amplicon ranges p13)
         (tssv)
         --minimum                        : $params.minimum # Reads a sequence needs
         --num_threads                    : $params.num_threads # Threads of FDSTOOLS tssv, at most the cores of the computer
         --baseQ                          : $params.baseQ # Minimum base quality: FDSTOOLS distrusts lower bases (0 = off)

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
    p00_parameters(log_text)
 
    // ──────────────── REFERENCE INDEXING ────────────────────
    p01_index_reference(params.reference)
    reference_index_ch = p01_index_reference.out.reference

    // ──────── HUMANS REFERENCE INDEXING (runs once) ─────────
    p02_index_humans(humans_index_ch, humans_base_ch)
    humans_index_ready_ch = p02_index_humans.out.humans_index_ready

    // ─────────── RAW READ MAPPING & FILTERING ───────────────
    p03_map_read_pairs(read_pairs_ch, params.reference, reference_index_ch)
    p04_remove_soft_clips(p03_map_read_pairs.out.raw_sam_ch, params.python_script_remove_scb, params.reference, reference_index_ch, params.left_primers, params.right_primers_rc)

    // ──────────────── MERGING & TRIMMING ────────────────────
    p05_merge_read_pairs(p04_remove_soft_clips.out.fastq_wo_scb_ch)
    p06_trim_primers(p05_merge_read_pairs.out, params.left_primers, params.right_primers_rc)

    // ────────────────── FINAL MAPPING ───────────────────────
    p07_map_merged_reads(p05_merge_read_pairs.out, params.reference, reference_index_ch, params.amplicon_middle_positions)
    p08_map_trimmed_reads(p06_trim_primers.out, params.reference, reference_index_ch, params.amplicon_middle_positions)

    // ────────────────── NUMTs FILTERING ─────────────────────
    // rtn runs once, on the primer-trimmed reads (p09); p10 keeps the same reads for FDSTOOLS
    // rtn takes the humans folder itself and the file name from p02 once the index is ready
    p09_filter_numts(p08_map_trimmed_reads.out, params.amplicon_middle_positions,  humans_index_ch, humans_index_ready_ch, numts_index_ch, numts_base_ch,
                                            rtn_lowq, params.python_script_rtn_qualities)
    p10_reads_for_fdstools(p07_map_merged_reads.out.join(p09_filter_numts.out.kept_names, by: 0))
    
    // ────────────────── QUALITY CONTROL ─────────────────────
    p11_read_depth(p09_filter_numts.out.bams, params.python_script_generate_read_depth_plot)
    

    // ────────────────── VARIANT CALLING ─────────────────────
    fdstools_ch = p12_fdstools(p10_reads_for_fdstools.out, params.fdstools_library)


    // ────────────── PROCESS & MERGE VARIANTS ────────────────
    p13_report(fdstools_ch, params.reference, params.python_script_process_fdstools_sast, params.python_script_merge_fdstools_mutect2,
                               files("$baseDir/resources/scripts/{notation,frames,iupac}.py"))

}
