# LD map

We compute some files representing linkage disequilibrium (LD) structure. 
Specifically, we generate the following for each population for each dataset type.

- `bool.prune.in`
- `bool.orune.out`
- `ld_map.ld.gz`

The first two files are the results of LD pruning. 

The last one is the output from `--r2`.

## reformatting output from `r2`

The output from r2 is a table file, but with unclear column separater.

```{bash}
cd /oak/stanford/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/ldmap
zcat ukb24983_ukb24983_cal_hla_imp.white_british.ld_map.ld.gz | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' | sed -e "s/^CHR_A/#CHR_A/g" | bgzip -l9 -@6 > ukb24983_ukb24983_cal_hla_imp.white_british.ld_map.tsv.gz
tabix -c '#' -s 1 -b 2 -e 5 ukb24983_ukb24983_cal_hla_imp.white_british.ld_map.tsv.gz
```

We create tabix index for `CHR_A:[SNP_A, BP_B]` interval.

## datasets

- `array_imp_combined_no_cnv`: the LD map for the combined dataset of array (cal) + HLA + imputation (v3) (no CNV).
  - The results files are copied to the Google Drive folder: https://drive.google.com/drive/folders/1mwKZtGfrOadASYVui6FQcy8k-Kk7f5k5
- `imp_v3_ldmap`: an initial attempt to make LD map for the imputed dataset. This is not useful.

## `LD_lookup.sh`

By default, we will look up the LD file computed for `white_british` cohort using `array_imp_combined_no_cnv` dataset.

```{bash}
$ bash LD_lookup.sh --r2 0.95 1 838732
#CHR_A  BP_A    SNP_A   CHR_B   BP_B    SNP_B   R2
1       838732  1:838732:G:A    1       838665  1:838665:T:C    0.996646
1       838732  1:838732:G:A    1       838890  1:838890:A:G    0.996646
1       838732  1:838732:G:A    1       838916  1:838916:A:T    0.988727
1       838732  1:838732:G:A    1       839461  1:839461:T:C    0.997124
1       838732  1:838732:G:A    1       839495  1:839495:G:T    0.999588
1       838732  1:838732:G:A    1       839528  1:839528:A:G    0.994058
1       838732  1:838732:G:A    1       839529  1:839529:T:G    0.996713
```
