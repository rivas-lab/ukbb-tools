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

## check the completion of the GWAS

https://github.com/rivas-lab/ukbb-tools/issues/19

```
cat 4_combined.idx.tsv | awk '(NR>1){print $3}' | while read f ; do echo $(zcat $f | wc -l) $f ; done | tr ' ' '\t' | tee /dev/stderr > 5_counts.$(date +%Y%m%d-%H%M%S).tsv
```

```
$ cat 5_counts.20200711-131359.tsv | awk '$1 != 1080278'
1080969 /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/white_british/ukb24983_v2_hg19.INI22152.array-combined.glm.linear.gz
```

The line counts indicates the file is now fixed.

This patch generated 824 files with 1080278 lines (-691) because the chrY variants were skipped and one file with 1080969 lines.
