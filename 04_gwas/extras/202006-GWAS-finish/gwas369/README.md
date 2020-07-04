# gwas 369

`/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/2007183/40831/others/ukb24983_v2_hg19.INI22149.array-combined.glm.linear.gz`

This file has 1080600 lines, meaning that there are 369 variants missing.

We first generate the list of missing variants.

```{bash}
1_list_missing_vars.sh
```

We then check the overlap with the `one_array` list.

```{bash}
bash 2_one_array.sh missing.369.lst
```

It turned out that all variants are in both arrays.

We then run GWAS.

```{bash}
bash 3_gwas.sh
```

We then combine the results

```{bash}
bash 4_combine.sh /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/2007183/40831/others/ukb24983_v2_hg19.INI22149.array-combined.glm.linear.gz  data/others/ukb24983_v2_hg19.INI22149.array-combined.PHENO1.glm.linear
```

The line counts indicates the file is now fixed.
