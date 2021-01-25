# LD map and LD indep flag for the combined dataset of array and exome 200k

## Yosuke Tanigawa, 2021/1/12

We previously prepared [array & exome combined dataset (hg19)](/03_filtering/array-exome-combined/20201217/) and here we compute the LD map.

For 4 populations (`non_british_white`, `african`, `s_asian`, and `e_asian`), this is a very straightforward task. We have [`1_ldmap.sbatch.sh`](1_ldmap.sbatch.sh).

For `white_british`, we split the task for each chromosome and run `plink1.9 --r2` for each chromosome separately to reduce the computational time (wall time) and concatenated the resulting files. We use [`1_ldmap-WB.sbatch.sh`](1_ldmap-WB.sbatch.sh) and [`1_ldmap-WB-combine.sh`](1_ldmap-WB-combine.sh) for this analysis.

We subsequently performed the LD independence analysis with [the procedure similar to what we used in the array-combined dataset](/17_annotation/20201012_array-combined). Specifically, we used an iterative procedure to apply `plink2 --indep-pairwise 1000kb 1 0.5` and prioritize the variants in the following order:

1. The [LD pruned set of variants in the array combined dataset](/17_annotation/20201012_array-combined)
  - This itself prioritized PTVs, PAVs, PCVs, UTR-region variants, and intronic variants, over the remaining variants in this order.
2. We removed all the variants that are in LD with the previously selected set.
3. We performed the LD pruning for the remaining (exome 200k) variants with priority assigned to PTVs, PAVs, PCVs, UTR-region variants, intronic variants, and the remaining variants, in this order.

We incorporated the results of the LD pruning into [the variant annotation file](/03_filtering/array-exome-combined/20210111).


## versino history

### `20210112`

- The proper version of the LD map for the array-exome-combined dataset
- The results were integrated into [the variant annotation file](/03_filtering/array-exome-combined/20210111), `ukb24983_cal_hla_cnv_exomeOQFE.annot_20210111.tsv`.

### `20201223` (please don't use this version)

- This is the initial version of the LD map computation
- It turned out that there was an error in the merged genotype dataset used in `ldmap_20201223` version.
  - Specifically, hg19/hg38 were not properly handled, which had a critical influence on the distance-based window (`--ld-window-kb 1000`).
- We revised the analysis.

## file locatoin

The results are stored in `/oak/stanford/groups/mrivas/ukbb24983/array-exome-combined/ldmap/ldmap_20210112`.

- The results from LD pruning procedure:
  - `ukb24983_cal_hla_exomeOQFE.<population>.bool.prune.in`
  - `ukb24983_cal_hla_exomeOQFE.<population>.bool.prune.out`
  - `ukb24983_cal_hla_exomeOQFE.<population>.bool.log`: plink log file
- The results from `plink --r2 gz --ld-window-r2 0.1 --ld-window-kb 1000 --ld-window 100000000` command:
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.tsv.gz`: the 7 column file containing the r2 value
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.tsv.gz.tbi`: the tabix index
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.log`: plink log file
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.hh`
  - `ukb24983_cal_hla_exomeOQFE.<population>.ld_map.nosex`

## Notes on analysis scripts

```{bash}
for csq in ptv pav pcv utr intron others ; do echo "[$(date +%Y%m%d-%H%M%S)] $csq" ; bash 2c_LD_indep.sh $csq ;  done
```

