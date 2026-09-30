#!/bin/bash
# Run a haplogroup list in batches with run_batch.sh. Finished batches are
# skipped, so the run can be restarted; stops when free disk space drops
# below MIN_FREE_GB (default 20).
#
# Usage: MOTIFS=<mitoLEAF hgmotifs.json> run_all.sh <haplogroup_list> <out_dir> [batch_size] [reads_per_amplicon]
set -euo pipefail

LIST=$1
OUT=$(mkdir -p "$2" && cd "$2" && pwd)
SIZE=${3:-100}
READS=${4:-100}
MIN_FREE_GB=${MIN_FREE_GB:-20}
HERE=$(cd "$(dirname "$0")" && pwd)

mkdir -p "$OUT/lists"
[ -n "$(ls "$OUT/lists")" ] || split -l "$SIZE" -a 3 "$LIST" "$OUT/lists/batch_"

for list in "$OUT"/lists/batch_*; do
    batch=$(basename "$list")
    [ -f "$OUT/$batch/keep/score/summary.tsv" ] && continue
    free=$(df -Pk "$OUT" | awk 'NR==2 {print int($4 / 1048576)}')
    if [ "$free" -lt "$MIN_FREE_GB" ]; then
        echo "$(date '+%F %H:%M') stopped: only ${free} GB free" >> "$OUT/progress.txt"
        exit 1
    fi
    if "$HERE/run_batch.sh" "$list" "$OUT/$batch" "$READS"; then
        echo "$(date '+%F %H:%M') $batch done" >> "$OUT/progress.txt"
    else
        echo "$(date '+%F %H:%M') $batch FAILED (work folder kept)" >> "$OUT/progress.txt"
    fi
done
