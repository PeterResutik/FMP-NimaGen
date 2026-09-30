# Simulated haplogroups

Tools to run the pipeline on simulated NimaGen reads for mitoLEAF haplogroups
and compare the calls with the known truth.

Haplogroup motifs come from mitoLEAF (`docs/data/hgmotifs.json` in
https://github.com/forensicgenomics/mitoLeaf, MPL-2.0), which is not copied
into this repository.

- `simulate_reads.py` builds each haplogroup's genome (rCRS plus its motif;
  the rCRS `N` at 3107 is absent) and writes read pairs as in real MiSeq
  data: 226 bases from each end (`--read-length`), primers from the primer
  files, inserts located by the primers on rCRS, substitution errors from
  per-base qualities. Positions mitoLEAF marks as ambiguous (IUPAC codes
  such as `16519Y`, lowercase entries such as `315.1c`) get one of their
  states, drawn per haplogroup (`--ambiguous one`, the default), become
  50/50 heteroplasmies on two haplotypes (`--ambiguous het`), or always take
  the variant (`--ambiguous alt`) or the rCRS state (`--ambiguous rcrs`).
  Drawn states and `alt` can combine in ways few real genomes show (e.g.
  `16181C 16182C 16183C` without `16189C`), which rtn then removes. No PCR
  stutter yet.
- `score_calls.py` compares the merged Excel with the truth, overall and per
  caller, and lists mismatches (C-stretches 303-315 and 16180-16193 separately).
  `--exclude` leaves out windows, e.g. `300-320,16170-16200` or `16519-16519`.
- `run_batch.sh` simulates a list of haplogroups, runs the pipeline with
  `--skip_numt_filter --publish_bams false` (`NUMT_FILTER=1` keeps rtn),
  keeps the small outputs
  (`sast.csv`, Mutect2 VCF, merged Excel, truth, parameter log), scores them
  and deletes reads, work folder and full output. `QUEUE_SIZE` (default 6)
  and `SIM_PROCESSES` (default 8) set the parallelism.
- `run_all.sh` runs a long list in batches, skips finished batches and stops
  when disk space runs low.

```
conda activate FMP-NimaGen
MOTIFS=path/to/mitoLeaf/docs/data/hgmotifs.json \
  tests/simulation/run_batch.sh haplogroups.txt /path/to/batch_dir 100
```
