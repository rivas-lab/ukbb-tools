# Exome 200k GWAS w/ array job script

## Yosuke Tanigawa, 2020/10/26

One need to be very careful when applying the script across multiple computing cluster systems.

## scripts

- [`0_functions.sh`](0_functions.sh): functions
- [`1_plink.gwas.sh`](1_plink.gwas.sh): run GWAS for a subset of variants

The rest is still under development

- [`2_check_results.sh`](2_check_results.sh): check if the output files
- [`3_combine_results.sh`](3_combine_results.sh): combine the output files

## usage (copied from the GWAS script for the array-combined script)

By default, the script apply GWAS for the array dataset and split the task into 100 chunks.
You can run each run with a specified index (`${idx}`) for `idx=1..100`, check the results, and combine the files.

```{bash}
bash 1_plink.gwas.sh      --GBE_ID ${GBE_ID} ${idx} ${output_dir}
bash 2_check_results.sh   --GBE_ID ${GBE_ID} ${output_dir}
bash 3_combine_results.sh --GBE_ID ${GBE_ID} ${output_dir}
```


