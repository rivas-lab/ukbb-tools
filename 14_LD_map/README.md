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

```
cd /oak/stanford/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/ldmap
zcat ukb24983_ukb24983_cal_hla_imp.white_british.ld_map.ld.gz | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' | sed -e "s/^CHR_A/#CHR_A/g" | bgzip -l9 -@6 > ukb24983_ukb24983_cal_hla_imp.white_british.ld_map.tsv.gz
tabix -c '#' -s 1 -b 2 -e 5 ukb24983_ukb24983_cal_hla_imp.white_british.ld_map.tsv.gz
```

We create tabix index for `CHR_A:[SNP_A, BP_B]` interval.

## datasets

- `array_imp_combined_no_cnv`: the LD map for the combined dataset of array (cal) + HLA + imputation (v3) (no CNV).
  - The results files are copied to the Google Drive folder: https://drive.google.com/drive/folders/1mwKZtGfrOadASYVui6FQcy8k-Kk7f5k5
- `imp_v3_ldmap`: an initial attempt to make LD map for the imputed dataset. This is not useful.

