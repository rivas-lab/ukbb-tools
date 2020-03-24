# convert the BGEN files to pgen files (imputation v3 dataset)

Those scripts are originally written and executed in 
`/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/`

## TR;DL How to regenerate pgen files in SCRATCH?

```
sbatch --array=1-24 bgen2pgen.sbatch.sh
```


## Step 1. conversion of BGEN to pgen

`bgen2pgen.sbatch.sh`

## Step 2. variant ID

`plink_var_ids`

We converted the file format and assigned variant IDs based on CHR:POS:REF:ALT

We performed a further manual curation for `PAR`.

```
[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/ukbb24983/imp/pgen]$ zstdcat ../pgen/ukb24983_imp_chrXY_v3.pvar.zst | grep -v '#' | cut -f 1| uniq -c
      1 PAR1
  45905 XY

[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/ukbb24983/imp/pgen]$ mv ukb24983_imp_chrXY_v3.pvar.zst ukb24983_imp_chrXY_v3.tmp.pvar.zst ; zstdcat ukb24983_imp_chrXY_v3.tmp.pvar.zst | grep -A1 PAR | sed -e "s/PAR1/XY/g" | zstd -15 > ukb24983_imp_chrXY_v3.pvar.zst

[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/ukbb24983/imp/pgen]$ zstdcat ../pgen/ukb24983_imp_chrXY_v3.tmp.pvar.zst | head
##INFO=<ID=ORIGINAL_VAR_ID,Number=1,Type=String,Description="The variant ID in the original data release">
#CHROM  POS     ID      REF     ALT     INFO
PAR1    60014   PAR1:60014:T:C  T       C       ORIGINAL_VAR_ID=rs370048753
XY      60014   XY:60014:T:G    T       G       ORIGINAL_VAR_ID=rs370048753
XY      60017   XY:60017:C:T    C       T       ORIGINAL_VAR_ID=X:60017_C_T
XY      60060   XY:60060:G:C    G       C       ORIGINAL_VAR_ID=rs148832940
XY      60072   XY:60072:G:C    G       C       ORIGINAL_VAR_ID=rs116895855
XY      60072   XY:60072:G:T    G       T       ORIGINAL_VAR_ID=rs116895855
XY      60112   XY:60112:G:C    G       C       ORIGINAL_VAR_ID=rs111065979
XY      60146   XY:60146:G:C    G       C       ORIGINAL_VAR_ID=rs138058540

[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/ukbb24983/imp/pgen]$ zstdcat ../pgen/ukb24983_imp_chrXY_v3.pvar.zst | head
##INFO=<ID=ORIGINAL_VAR_ID,Number=1,Type=String,Description="The variant ID in the original data release">
#CHROM  POS     ID      REF     ALT     INFO
XY      60014   XY:60014:T:C    T       C       ORIGINAL_VAR_ID=rs370048753
XY      60014   XY:60014:T:G    T       G       ORIGINAL_VAR_ID=rs370048753
XY      60017   XY:60017:C:T    C       T       ORIGINAL_VAR_ID=X:60017_C_T
XY      60060   XY:60060:G:C    G       C       ORIGINAL_VAR_ID=rs148832940
XY      60072   XY:60072:G:C    G       C       ORIGINAL_VAR_ID=rs116895855
XY      60072   XY:60072:G:T    G       T       ORIGINAL_VAR_ID=rs116895855
XY      60112   XY:60112:G:C    G       C       ORIGINAL_VAR_ID=rs111065979
XY      60146   XY:60146:G:C    G       C       ORIGINAL_VAR_ID=rs138058540
```

The original variant IDs are kept in the INFO fields

## Step 3. Concatenate the information in the `mfi` files

`concatenate_mfi.sh`

The resulting file has allele frequency (A1) and INFO score.

