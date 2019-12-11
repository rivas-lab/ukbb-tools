# The genetic relationship matrix (GRM)

Using GCTA, we compute the genetic relationship matrix (GRM) using the directly genotyped data in the UK Biobank (data release version 2).

https://cnsgenomics.com/software/gcta/#Overview

## Data availability

The results files are in `/oak/stanford/groups/mrivas/ukbb24983/cal/grm`

- `ukb24983_cal_cALL_v2_hg19.white_british.grm.bin`
- `ukb24983_cal_cALL_v2_hg19.white_british.grm.id`
- `ukb24983_cal_cALL_v2_hg19.white_british.grm.N.bin`
- `ukb24983_cal_cALL_v2_hg19.white_british.log`

## Input data

We used the directly genotyped data (`/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19*`).
For computatinoal efficiencly, we used a copy in `/scratch`.

## Methods

We used `--make-grm-part` command to compute the GRM matrix using 100 jobs and combied them together.

- `15_GRM_misc.sh`: this script defines the GCTA commands used in the computation.
- `gcta_grm.part.sbatch.sh`: This SBATCH script was used to run `--make-grm-part` commands. (`$ sbatch --array=1-100 gcta_grm.part.sbatch.sh`)
- `gcta_grm.part.combine.sh`: This script combines the results into one file.
