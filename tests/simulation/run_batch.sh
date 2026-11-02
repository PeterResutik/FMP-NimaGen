#!/bin/bash
# Simulate a list of haplogroups, run the pipeline on them, keep only the small
# outputs (sast.csv, merged Excel, truth, parameter log), score
# them against the truth and delete reads, work folder and full output.
# Before the work folder goes, the report step is run again at each detection
# threshold in EXTRA_THRESHOLDS (default "3 4 5"; empty for none): only p13
# reruns, and keep/merged_<t> and keep/score_<t> hold its reports and scores
# (keep/merged and keep/score are those of the pipeline's defaults).
#
# Usage: MOTIFS=<mitoLEAF hgmotifs.json> run_batch.sh <haplogroup_list> <batch_dir> [reads_per_amplicon]
# With MIXTURES=1 the list holds mixtures ("major minor share" per line,
# simulate_reads.py --mixtures), scored by score_mixtures.py.
# Needs the FMP-NimaGen conda environment to be active. rtn is skipped
# (simulated reads contain no NUMTs) unless NUMT_FILTER=1. QUEUE_SIZE (default 6)
# sets how many pipeline tasks run at once. ERROR_PROFILE (made by
# build_error_profile.py) gives the reads the qualities and errors of real
# reads. On failure the work folder is kept for debugging.
set -euo pipefail

LIST=$1
BATCH=$(mkdir -p "$2" && cd "$2" && pwd)
READS=${3:-100}
MOTIFS=${MOTIFS:?set MOTIFS to mitoLEAF docs/data/hgmotifs.json}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=${REPO:-$(cd "$HERE/../.." && pwd)}

INPUT=(--haplogroup-list "$LIST")
[ "${MIXTURES:-0}" = 1 ] && INPUT=(--mixtures "$LIST")
python "$HERE/simulate_reads.py" --motifs "$MOTIFS" "${INPUT[@]}" \
    --outdir "$BATCH/reads" --reads-per-amplicon "$READS" --processes "${SIM_PROCESSES:-8}" \
    --reference "$REPO/resources/rCRS/rCRS_NimaGen.fasta" \
    --left-primers "$REPO/resources/primers/left_primers.fasta" \
    --right-primers "$REPO/resources/primers/right_primers_rc.fasta" \
    ${ERROR_PROFILE:+--error-profile "$ERROR_PROFILE"} > "$BATCH/simulate.log"

FLAGS=(--publish_bams false)
[ "${NUMT_FILTER:-0}" = 1 ] || FLAGS+=(--skip_numt_filter)

# a run name of our own, to resume exactly this run for the other thresholds
# (other batches may run from the same folder at the same time)
RUN="sim_$(basename "$BATCH" | tr -c 'a-zA-Z0-9_\n' _)_$(date +%s)_$$"
start=$(date +%s)
(cd "$REPO" && nextflow run main.nf -name "$RUN" -profile local -qs "${QUEUE_SIZE:-6}" -w "$BATCH/work" \
    --reads "$BATCH/reads/*_{R1,R2}_001.fastq.gz" --outdir "$BATCH/out" \
    "${FLAGS[@]}") > "$BATCH/nextflow.log" 2>&1
echo "pipeline seconds: $(( $(date +%s) - start ))" > "$BATCH/timing.txt"
SESSION=$(awk -F'\t' -v n="$RUN" '$3 == n {print $6}' "$REPO/.nextflow/history" | tail -1)
[ -n "$SESSION" ] || { echo "run $RUN not found in $REPO/.nextflow/history" >&2; exit 1; }

KEEP="$BATCH/keep"
mkdir -p "$KEEP/sast" "$KEEP/merged" "$KEEP/truth"
cp "$BATCH"/out/p12_fdstools/*/*.sast.csv "$KEEP/sast/"
cp "$BATCH"/out/p13_report/*.xlsx "$KEEP/merged/"
cp "$BATCH"/reads/*_truth.tsv "$KEEP/truth/"
cp "$BATCH/out/p00_parameters.txt" "$KEEP/"

score() { # merged_dir out_dir
    local scorer=score_calls.py
    [ "${MIXTURES:-0}" = 1 ] && scorer=score_mixtures.py
    python "$HERE/$scorer" --truth-dir "$KEEP/truth" --merged-dir "$1" --out "$2" \
        --reference "$REPO/resources/rCRS/rCRS_NimaGen.fasta"
}
score "$KEEP/merged" "$KEEP/score"

for t in ${EXTRA_THRESHOLDS-3 4 5}; do
    start=$(date +%s)
    (cd "$REPO" && nextflow run main.nf -profile local -resume "$SESSION" -qs "${QUEUE_SIZE:-6}" -w "$BATCH/work" \
        --reads "$BATCH/reads/*_{R1,R2}_001.fastq.gz" --outdir "$BATCH/out_$t" \
        --min_vf_FDS "$t" "${FLAGS[@]}") > "$BATCH/nextflow_$t.log" 2>&1
    echo "report step at $t%, seconds: $(( $(date +%s) - start ))" >> "$BATCH/timing.txt"
    mkdir -p "$KEEP/merged_$t"
    cp "$BATCH/out_$t"/p13_report/*.xlsx "$KEEP/merged_$t/"
    score "$KEEP/merged_$t" "$KEEP/score_$t"
done

rm -rf "${BATCH:?}/work" "${BATCH:?}/reads" "${BATCH:?}/out" "${BATCH:?}"/out_*
