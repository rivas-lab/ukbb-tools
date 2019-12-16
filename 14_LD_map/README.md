# LD map

We compute some files representing linkage disequilibrium (LD) structure. 
Specifically, we generate the following for each population for each dataset type.

- `bool.prune.in`
- `bool.orune.out`
- `ld_map.ld.gz`

The first two files are the results of LD pruning. 

The last one is the output from `--r2`.

## datasets

- `array_imp_combined_no_cnv`: the LD map for the combined dataset of array (cal) + HLA + imputation (v3) (no CNV).
  - The results files are copied to the Google Drive folder: https://drive.google.com/drive/folders/1mwKZtGfrOadASYVui6FQcy8k-Kk7f5k5
- `imp_v3_ldmap`: an initial attempt to make LD map for the imputed dataset. This is not useful.

