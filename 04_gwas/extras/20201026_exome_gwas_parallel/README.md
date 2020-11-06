# Exome 200k GWAS scan

## Yosuke Tanigawa, 2020/10/26-2020/11/5

We apply GWAS scan for Exome 200k dataset.

## ToDo (as of 2020/11/5)

- copy/move the files from `/scratch` to `/oak`
- meta-analysis. we currently apply `ERRCODE=='.'` filter.

## data location

We tentatively have the data in `/scratch`, `/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE`.

## GWAS analysis

We used master phe file (version `20201002`) for quantitative phenotypes (`INI` and `QT_FC`) and phe files pointed in [`phenotype_info.tsv`](/05_gbe/array-combined/phenotype_info.tsv) for the binary phenotypes.
We analyzed 7 population groups (`white_british`, `non_british_white`, `african`, `s_asian`, `e_asian`, `related`, and `others`) based on the keep file version 20200828.

To maximize the utilization of cluster computing resources, we split the genotype file (of 17M variants) into smaller pieces, apply GWAS for each chunk, and combine the results.
Specifically, we split the 17M variants into equally-sized 100 pieces for most analysis, with the exception of GWAS analysis of quantitative traits for white British individuals.

We used two computing clusters, Sherlock and scg. Because `/scratch` is not mounted on `scg`, we dumped some results under `oak` space. We moved/copied the results file from `oak` to `scratch` space so that we don't hit the disk uasge limits.

## Basic QC

We counted the number of lines in each file.

## Meta-analysis with Metal

We performed meta-analysis with `metal`.

## analysis scripts

We have several key analysis scripts and many wrapper scripts.

- [`0_functions.sh`](0_functions.sh): functions
- [`1_plink.gwas.sh`](1_plink.gwas.sh): run GWAS for a subset of variants
- [`4_combine_output.sh`](4_combine_output.sh): a script to merge files
- [`6a_metal.sh`](6a_metal.sh): meta-analysis script
- [`7a_wc_l.sh`](7a_wc_l.sh): a script to run `wc -l` for basic QC.
- [`8_patch_fix.sh`](8_patch_fix.sh): a helper script to fix issues in the intermediate (truncated) GWAS results

### Some example usage of the analysis scripts

#### merge the GWAS results in chunks

```{bash}
bash 3a_check_results.BINs.sh
bash 3b_merge_job_list.sh 3a_check_results.BINs.20201029-185932.tsv
wc 3b_merge_job_list.20201029-185957.tsv
sbatch -p mrivas --qos=high_p --time=1:0:00 --mem=6000 --nodes=1 --cores=1 --job-name=combine --output=logs_scratch/combine.%A_%a.out --error=logs_scratch/combine.%A_%a.err --array=1-57 ${parallel_sbatch_sh} 4_combine_output_BIN.wrapper.sh ${parallel_idx} 5 3b_merge_job_list.20201029-185957.tsv
```

