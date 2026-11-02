#!/usr/bin/env bash
# Builds humans_NimaGen_2026-10-30.fa, the human mitogenomes rtn compares reads with, from public
# sources (see README.md). Needs bwa, samtools and the pipeline's conda environment, and downloads
# from GitHub (mitoLEAF) and NCBI (about 70,000 genomes, one to two hours).
#
# usage: build.sh workdir
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
RES=$(cd "$HERE/../../.." && pwd)
W=$1
mkdir -p "$W"
NAME=humans_NimaGen_2026-10-30.fa
NUMTS=$RES/rtn_files/numts/Calabrese_Dayama_Smart_Numts_modified.fa
MITOLEAF=https://raw.githubusercontent.com/forensicgenomics/mitoLeaf/25b6076/docs/data

echo "1. mitoLEAF tables (version 1.6.1, commit 25b6076)"
for f in profiles.csv:f7f00e93575fc70e7052330dfc259437 hgmotifs.json:39521b01741850ece285b4806715bbb7 \
         tree.json:402b206e62759e2671cdf91ca86af85a; do
    [ -s "$W/${f%%:*}" ] || curl -sfL "$MITOLEAF/${f%%:*}" -o "$W/${f%%:*}"
    python -c "import hashlib, sys; sys.exit(hashlib.md5(open(sys.argv[1], 'rb').read()).hexdigest() != sys.argv[2])" \
        "$W/${f%%:*}" "${f#*:}" || { echo "${f%%:*}: md5 differs from the version used"; exit 1; }
done

echo "2. GenBank genomes on mitoLEAF's list"
python -c "import csv, sys; [print(r['accession']) for r in csv.DictReader(open(sys.argv[1])) if r['source'] == 'NCBI']" \
    "$W/profiles.csv" > "$W/accessions.txt"
python "$HERE/fetch_genbank.py" "$W/accessions.txt" "$W/genbank"
cat "$W"/genbank/batch_*.fa > "$W/genbank.fa"

echo "3. rCRS without the N at 3107"
python -c "import sys; s = ''.join(l.strip() for l in open(sys.argv[1]) if l[0] != '>')[:16569]; print('>chrM'); print(s[:3106] + s[3107:])" \
    "$RES/rCRS/rCRS_NimaGen.fasta" > "$W/ref.fa"

echo "4. screen the GenBank genomes"
mkdir -p "$W/scan" && cp "$W/ref.fa" "$W/scan/"
python "$HERE/genbank_singles.py" "$W/genbank.fa" "$W/scan"
python "$HERE/scan_humans.py" "$W/genbank.fa" "$W/hgmotifs.json" "$W/scan"
python "$HERE/amplicon_check.py" "$W/scan" "$NUMTS" "$RES"

echo "5. screen the re-called 1000 Genomes and HGDP people"
mkdir -p "$W/scan_people" && cp "$W/ref.fa" "$W/scan_people/"
python "$HERE/people_sequences.py" "$HERE/people_calls.tsv.gz" "$W/ref.fa" > "$W/scan_people/singles.fa"
python "$HERE/scan_humans.py" "$W/scan_people/singles.fa" "$W/hgmotifs.json" "$W/scan_people"
python "$HERE/amplicon_check.py" "$W/scan_people" "$NUMTS" "$RES"

echo "6. mitoLEAF haplogroup sequences from the screened genomes, and assemble"
python "$HERE/build_humans.py" "$W" "$HERE/people_calls.tsv.gz" "$W/$NAME"
bzip2 -kf "$W/$NAME"
md5sum "$W/$NAME" 2> /dev/null || md5 "$W/$NAME"
echo "done: $W/$NAME.bz2 (copy to ../$NAME.bz2), left_out.tsv, members.tsv and haplogroup_marks.tsv (copy here)"
