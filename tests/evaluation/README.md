# Evaluation on real data

Scripts that evaluate the pipeline on real samples. The simulations are in
`../simulation/`. No data are included: the scripts take the pipeline's
output folders and the reads as arguments.

- `mix_reads.py` mixes the read pairs of two real samples at a known share,
  amplicon by amplicon: every amplicon gets exactly the given share of its
  read pairs from the second sample, whatever the two samples' depths per
  amplicon. A pair's amplicon is the left primer its R1 starts with; pairs
  without one (primer dimers) form one more group. Read names get the prefix
  `mixA_` or `mixB_`, so each sample's reads can be counted after the
  pipeline. Name the output like `mix-<A>-<B>-m<share>_L001`.
- `evaluate_mixtures.py` compares the mixtures' reports with what the two
  samples' own reports predict. At every position where their major bases
  differ (except where either sample has a heteroplasmy, insertion or
  deletion of its own, and outside the C-stretches), the expected share
  comes from the `mixA_`/`mixB_` reads in
  the pipeline's BAM there; a share that reaches `--min-vf` is expected as a
  heteroplasmy (IUPAC code), a lower one as the major base. The mixtures need
  `--publish_bams true`, and the two samples need their own run with the
  same pipeline version.

```
python tests/evaluation/mix_reads.py A_R1_001.fastq.gz B_R1_001.fastq.gz 5 reads/mix-A-B-m5_L001
nextflow run main.nf --reads "reads/*_{R1,R2}_001.fastq.gz" --outdir out_mix --publish_bams true
python tests/evaluation/evaluate_mixtures.py out_mix out_pure --out eval_mix --min-vf 2
```
