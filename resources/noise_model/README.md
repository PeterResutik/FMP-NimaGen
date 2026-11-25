# Noise model

The depth rule (README, Depth Rule; `resources/scripts/depth_rule.py`) comes
from a noise model: for each control sample, how many wrong bases reach each
share of its reads, every 0.05% from 0.1 to 2% (with at least 2 reads). Each
line of the file is one control.

`2800M_8controls.tsv` holds the Promega 2800M positive controls of 8 MiSeq
runs, named by flow cell, run through version 0.6.0 of the pipeline with
`--minimum 1`. Left out: insertions and deletions, the C-stretches 57–60,
300–315 and 16180–16193, and 2800M's own low-level heteroplasmies at 1282,
2846, 8078, 10000, 12630, 14847 and 15759.

## A model from your own controls

1. Run the pipeline on the controls with `--minimum 1`, so that sequences seen
   once are kept:

   ```bash
   nextflow run main.nf -profile local --reads 'controls/*_{R1,R2}_001.fastq.gz' \
       --outdir controls_min1 --minimum 1
   ```

2. Build the model from their FDSTOOLS sequences:

   ```bash
   python resources/noise_model/build_noise_model.py \
       --library resources/fdstools/mtNG_library_file.txt \
       --reference resources/rCRS/rCRS_NimaGen.fasta \
       controls_min1/p12_fdstools/*/*.sast.csv > my_controls.tsv
   ```

   In each amplicon the most common sequence counts as the control's own, and
   every other substitution as a wrong base. A control DNA other than 2800M
   needs `--mask` with its own heteroplasmies (positions, comma-separated). A
   share of 0.1% needs 2 reads in 2,000, so a control with fewer reads per
   amplicon undercounts the lowest shares; the script says so. Use controls of
   several runs.

3. Give the model to the pipeline, alone or with the shipped controls; the
   controls of several files are pooled, each counting once:

   ```bash
   --noise_model my_controls.tsv
   --noise_model resources/noise_model/2800M_8controls.tsv,my_controls.tsv
   ```
