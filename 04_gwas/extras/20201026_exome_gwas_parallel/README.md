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

## Some script usage

```{bash}
find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/white_british-batch/logs/ -name "*.glm.log" | grep -v QT_ALL | parallel --eta -j6 -k 'bash 2c_mv_BINs_files_from_oak_to_scratch.sh {}'
bash 3a_check_results.BINs.sh
bash 3b_merge_job_list.sh 3a_check_results.BINs.20201029-185932.tsv
wc 3b_merge_job_list.20201029-185957.tsv
sbatch -p mrivas --qos=high_p --time=1:0:00 --mem=6000 --nodes=1 --cores=1 --job-name=combine --output=logs_scratch/combine.%A_%a.out --error=logs_scratch/combine.%A_%a.err --array=1-57 ${parallel_sbatch_sh} 4_combine_output_BIN.wrapper.sh ${parallel_idx} 5 3b_merge_job_list.20201029-185957.tsv

```

