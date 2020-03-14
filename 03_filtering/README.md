# Filtering

We do not employ sample- and variant-level QC in the traditional sense in this pipeline. Sample-level QC is often done by the BioBank itself before the data is sent out, and variant-level QC is done after GWAS is run (in case any signal can be recovered despite bad data).

The notebooks in this folder do two distinct tasks:

1) [Go over the overall marker quality of different BioBank arrays](https://github.com/rivas-lab/ukbb-tools/blob/master/03_filtering/Marker_QC.ipynb) (marker_QC.ipynb)
2) [Redact individuals from the dataset, define populations, and generate the GWAS covariate file](https://github.com/rivas-lab/ukbb-tools/blob/master/03_filtering/sample_qc_v3.1.ipynb) (rsample_qc_v3.1.ipynb)

## Sample level QC

### Participant Withdrawal

The latest set of individuals who withdrew their participation to UK Biobank is recorded in the following file:

`/oak/stanford/groups/mrivas/ukbb24983/sqc/W24983_20200204.csv`

This file is a superset of the previous versions of the withdrawal files.

#### Previous versions of withdrawal files

So far, we received the following 3 files.

- `w24983_20200204.csv` (117 individuals)
- `w24983_20181016.csv` (79 individuals)
- `w2498_20170726.csv` (3 individuals)

```{bash}
$ cd /oak/stanford/groups/mrivas/ukbb24983/sqc
$ comm -13 <( sort w24983_20200204.csv ) <( sort w24983_20181016.csv ) | wc -l
0
$ comm -13 <( sort w24983_20181016.csv ) <( sort w2498_20170726.csv ) | wc -l
0
```

### Population definition

#### summary

We used a combination of PCA (on array genotype data) and self-reported ancestry to define the following five population groups

- White British (`ukb24983_white_british.phe`, N = 337,138)
- Non-British White (`ukb24983_non_british_white.phe`, N = 24,905)
- South Asian (`ukb24983_s_asian.phe`, N = 7,885)
- African (`ukb24983_african.phe`, N = 6,497)
- East Asian (`ukb24983_e_asian.phe`, N = 1,154)

In total, we have 378,292 unrelated individuals

- Those files are available: `/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313`
  - In the directory, we have 

#### sample-level QC criteria

Our filtering criteria are as follows:

- Not present in the "remove file" (`! in_remove_file`, see 'Participant Withdrawal' section above)
- `FID >= 0`
- `IID >= 0`
- `putative_sex_chromosome_aneuploidy == 0`
- `het_missing_outliers == 0`
- `excess_relatives == 0`
- `used_in_pca_calculation == 1`

#### Genotype PC-based population definition

We defined the following thresholds to define ethnic groups

- White British
  - `-20 <= PC1 <= 40 && -25 <= PC2 <= 10` (Global PCs provided by UK Biobank)
  - `in_white_British_ancestry_subset == 1` (in sample QC file)
- Non-British White
  - `-20 <= PC1 <= 40 && -25 <= PC2 <= 10`
  - Based on self-reported ancestry (UKB field 21000), the individual is White, but is not White British.
- African
  - `260 <= PC1        &&   50 <= PC2` (Global PCs provided by UK Biobank)
  - Based on self-reported ancestry (UKB field 21000), the individual is White, but is not any of the followings: White, Asian, Mixed, and Other Ethnic Groups
- South Asian
  - `40 <= PC1 <= 120 && -170 <= PC2 <= -80` (Global PCs provided by UK Biobank)
  - Based on self-reported ancestry (UKB field 21000), the individual is White, but is not any of the followings: White, Black, Mixed, and Other Ethnic Groups
  - We subsequently applied thresholds on local PCs (the ones re-computed for the initially assigned population group)
  - `0.02 <= PC1 <= 0.03 && -0.05 <= PC2 <= 0.02`
- East Asian
  - `130 <= PC1 <= 170 &&         PC2 <= -230` (Global PCs provided by UK Biobank)
  - Based on self-reported ancestry (UKB field 21000), the individual is White, but is not any of the followings: White, Black, Mixed, and Other Ethnic Groups
  - We subsequently applied thresholds on local PCs (the ones re-computed for the initially assigned population group)
  - `-0.01 <= PC1 <= 0.02 && -0.02 <= PC2 <= 0`

See more details in notebook (v3.1)

#### Population-specific PCA

We applied PCA for 4 populations (all, but White British) and characterized population-specific PCs.

#### Additional clean-up using population specific PCs
