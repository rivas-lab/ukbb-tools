# LD map for the `cal + hla + imp` dataset


## File location

`/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined_no_cnv/`

For each population (`${pop}`), we have the following files:

```
ukb24983_ukb24983_cal_hla_imp.${pop}.bool.log
ukb24983_ukb24983_cal_hla_imp.${pop}.bool.prune.in
ukb24983_ukb24983_cal_hla_imp.${pop}.bool.prune.out
ukb24983_ukb24983_cal_hla_imp.${pop}.ld_map.hh
ukb24983_ukb24983_cal_hla_imp.${pop}.ld_map.ld.gz
ukb24983_ukb24983_cal_hla_imp.${pop}.ld_map.log
ukb24983_ukb24983_cal_hla_imp.${pop}.ld_map.nosex
```

where, 

- `*bool*` is the results from LD pruning (`plink2 --indep-pairwise 50 5 0.5`)
- `*ld_map*` is the results from `--ld-window-kb 1000 --ld-window-r2 0.1 --r2 gz`

### Back-up in Google Drive

- https://drive.google.com/drive/folders/1mwKZtGfrOadASYVui6FQcy8k-Kk7f5k5 

### Copy of the data in GBE

- `/opt/biobankengine/GlobalBioBankEngineRepo/gbe_browser/static/ldmap`
