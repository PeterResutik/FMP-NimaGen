#!/bin/bash
# Simulate a list of haplogroups, run the pipeline on them, keep only the small
# outputs (sast.csv, Mutect2 VCF, merged Excel, truth, parameter log), score
# them against the truth and delete reads, work folder and full output.
#
# Usage: MOTIFS=<mitoLEAF hgmotifs.json> run_batch.sh <haplogroup_list> <batch_dir> [reads_per_amplicon]
# Needs the FMP-NimaGen conda environment to be active. rtn is skipped
# (simulated reads contain no NUMTs) unless NUMT_FILTER=1. QUEUE_SIZE (default 6)
# sets how many pipeline tasks run at once. On failure the work folder is kept
# for debugging.
set -euo pipefail

LIST=$1
BATCH=$(mkdir -p "$2" && cd "$2" && pwd)
READS=${3:-100}
MOTIFS=${MOTIFS:?set MOTIFS to mitoLEAF docs/data/hgmotifs.json}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${REPO:-$(cd "$HERE/../.." && pwd)}

python "$HERE/simulate_reads.py" --motifs "$MOTIFS" --haplogroup-list "$LIST" \
    --outdir "$BATCH/reads" --reads-per-amplicon "$READS" --processes "${SIM_PROCESSES:-8}" \
    --reference "$REPO/resources/rCRS/rCRS_NimaGen.fasta" \
    --left-primers "$REPO/resources/primers/left_primers.fasta" \
    --right-primers "$REPO/resources/primers/right_primers_rc.fasta" > "$BATCH/simulate.log"

FLAGS=(--publish_bams false)
[ "${NUMT_FILTER:-0}" = 1 ] || FLAGS+=(--skip_numt_filter)

start=$(date +%s)
(cd "$REPO" && nextflow run main.nf -profile local -qs "${QUEUE_SIZE:-6}" -w "$BATCH/work" \
    --reads "$BATCH/reads/*_{R1,R2}_001.fastq.gz" --outdir "$BATCH/out" \
    "${FLAGS[@]}") > "$BATCH/nextflow.log" 2>&1
echo "pipeline seconds: $(( $(date +%s) - start ))" > "$BATCH/timing.txt"

KEEP="$BATCH/keep"
mkdir -p "$KEEP/sast" "$KEEP/vcf" "$KEEP/merged" "$KEEP/truth"
cp "$BATCH"/out/p11_fdstools/*/*.sast.csv "$KEEP/sast/"
cp "$BATCH"/out/p12_mutect2/*/*.vcf.gz "$KEEP/vcf/"
cp "$BATCH"/out/p13_merged_variants_xlsx/*.xlsx "$KEEP/merged/"
cp "$BATCH"/reads/*_truth.tsv "$KEEP/truth/"
cp "$BATCH/out/p00_parameters.txt" "$KEEP/"

python "$HERE/score_calls.py" --truth-dir "$KEEP/truth" --merged-dir "$KEEP/merged" --out "$KEEP/score" \
    --reference "$REPO/resources/rCRS/rCRS_NimaGen.fasta"

rm -rf "${BATCH:?}/work" "${BATCH:?}/reads" "${BATCH:?}/out"
