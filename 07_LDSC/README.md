# 07 LD score regression (LDSC)

We compute genetic parameters by [LD score regression (LDSC)](https://github.com/bulik/ldsc) with GWAS summary statistics as input.

## LDSC helper scripts

We have LDSC helper scripts. Please see [`helpers`](helpers) for more info.

## LDSC results files (version 2020)

We convert the GWAS summary statistics from the latest GWAS freeze. We keep track of that activity in [the issue ticket #26](https://github.com/rivas-lab/ukbb-tools/issues/26) for the 2020 refresh.

### Summary statistics in LDSC munge format

The converted summary statistics is stored in this directory:

`/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc`

Each population group has a sub-directory.


## LDSC results files (old ~ 2018 version?)

The results of LDSC are written in the following directories:

- `$OAK/projects/ukbb-tools-ldsc/h2`: Estimated heritability ($h^2$)
- `$OAK/projects/ukbb-tools-ldsc/munged`: "Munged" summary statistics. This directory contains summary statistics files in a format that can be used for LDSC (In fact we call [LDSC's `munge_sumstats.py`](https://github.com/bulik/ldsc/blob/master/munge_sumstats.py)) to generate these "munged" files.
- `/oak/stanford/groups/mrivas/projects/h2-estimation/private_output/ukb_ld_scores/TWINSUK`: LD score files. We use these LD score files for White British samples in the UK Biobank. These files are computed from other project ([Genetic parameter estimation project in the Rivas lab](https://github.com/rivas-lab/genetic-parameter-estimation)).

## Relevant info

- [LD score regression (LDSC)](https://github.com/bulik/ldsc)
- [Genetic parameter estimation project in the Rivas lab](https://github.com/rivas-lab/genetic-parameter-estimation). Some of the codes are follked from this repository. (commit ID: `c909bd998b5a73742433c1d2293a0928f82947ea`).
- We now have [Docker image for LDSC](https://github.com/yk-tanigawa/docker-LDSC)!

## Contents of this directory

- `helper`: This sub-directory contains helper scripts. Most of them can be used from command-line.
- `misc`: This sub-directory contains a list of variants that we use in LDSC.
- `jobs`: This sub-directory contains example job files to run LDSC on Sherlock cluster. The set of scripts in this directory is a folk from other repo ([Yosuke's array job template](https://github.com/yk-tanigawa/array-job-template)).
- `scripts`: This sub-directory contains scripts to run LDSC.
