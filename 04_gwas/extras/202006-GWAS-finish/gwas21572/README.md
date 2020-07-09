# gwas 21572

- There are 880 files with 1059397 lines (meaning that 21572 lines are missing in each file).
- this script checks whether those files have the same set of variants

```
bash ../gwas691/1_check.nvars.sh 1059397
```

We confirmed that those 880 files have the same set of 1,059,397 variants.

We generate the list of missing variants.

```
bash 1b_list_missing_vars.sh 1_check.nvars.out.lst
```

This writes the 21,572 missing variants into `missing.lst`.

We then check the overlap with the `one_array` list.

```{bash}
bash ../gwas369/2_one_array.sh missing.lst
```

It turned out that

- 21,033 variants are on both arrays
- 539 variants are on one array

We then run GWAS.

Test case:

```{bash}
bash 3_gwas.sh BIN20477 african
```

It turned out that chrY variants were skipped.

```
Note: Skipping chrY in --glm regression on phenotype 'PHENO1', since #
remaining samples <= # predictor columns.
```

Some might fail but we applied batch run.

```{bash}
bash 3_gwas.batch.sh 1059397
```

@@@@ We are still running GWAS as of 2020/7/4 1 pm

We then combine the results

```{bash}
bash 4_combine.sh /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/2007183/40831/others/ukb24983_v2_hg19.INI22149.array-combined.glm.linear.gz  data/others/ukb24983_v2_hg19.INI22149.array-combined.PHENO1.glm.linear
```

The line counts indicates the file is now fixed.
