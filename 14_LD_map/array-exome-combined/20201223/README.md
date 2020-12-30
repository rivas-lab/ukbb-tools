# LD map for the combined dataset of array and exome 200k

## Yosuke Tanigawa, 2020/12/23-12/30

We previously prepared [array & exome combined dataset (hg19)](/03_filtering/array-exome-combined/20201217/) and here we compute the LD map.

For 4 populations (`non_british_white`, `african`, `s_asian`, and `e_asian`), this is a very straightforward task. We have [`1_ldmap.sbatch.sh`](1_ldmap.sbatch.sh).

For `white_british`, we split the task for each chromosome and run `plink1.9 --r2` for each chromosome separately to reduce the computational time (wall time) and concatenated the resulting files. We use [`1_ldmap-WB.sbatch.sh`](1_ldmap-WB.sbatch.sh) and [`1_ldmap-WB-combine.sh`](1_ldmap-WB-combine.sh) for this analysis.

## file locatoin

The results are stored in `/oak/stanford/groups/mrivas/ukbb24983/array-exome-combined/ldmap/ldmap_20201223`.

- The results from `plink2 --indep-pairwise 50 5 0.5` command:
  - `ukb24983_cal_hla_exomeOQFE.<population>.bool.prune.in`
  - `ukb24983_cal_hla_exomeOQFE.<population>.bool.prune.out`
  - `ukb24983_cal_hla_exomeOQFE.<population>.bool.log`: plink log file
- The results from `plink --r2 gz --ld-window-r2 0.1 --ld-window-kb 1000 --ld-window 100000000` command:
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.tsv.gz`: the 7 column file containing the r2 value
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.tsv.gz.tbi`: the tabix index
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.log`: plink log file
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.hh`
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.nosex`
