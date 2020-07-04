# gwas 691

- There are 407 files with 1080278 lines (meaning that 691 lines are missing in each file).
- this script checks whether those files have the same set of variants

```
bash 1_check.nvars.sh 1080278
```

We confirmed that those 407 files have the same set of 1,080,278 variants.


We generate the list of missing variants.

```
bash 1b_list_missing_vars.sh 1_check.nvars.out.lst
```

This writes the 691 missing variants into `missing.lst`.

We then check the overlap with the `one_array` list.

```{bash}
bash ../gwas369/2_one_array.sh missing.lst
```

It turned out that

- 369 variants are on both arrays
- 322 variants are on one array

We then run GWAS.

Test case:

```{bash}
bash 3_gwas.sh BIN20445 others
```

It turned out that the glm was skipped.

```
grep 'Skipping case/control phenotype' data/others/ukb24983_v2_hg19.BIN20445.array-combined.log
--glm: Skipping case/control phenotype 'PHENO1' since all samples are controls.
--glm: Skipping case/control phenotype 'PHENO1' since all samples are controls.
```

Some might fail but we applied batch run.

```{bash}
bash 3_gwas.batch.sh 1080278
```

We check the results:

```{bash}
cd /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas691

find data -name "*.log" | wc
    407     407   23387

find data -name "*glm*" | wc
    180     180   13457
```

This indicates the computation was applied to 407 (pop, GBE_ID) pairs and obtained the results for 180 of them.

```
$ cat 4_wc_l.tsv | awk '$1 != 692'
#wc_l   file
323     data/others/ukb24983_v2_hg19.INI22154.array-combined.PHENO1.glm.linear
323     data/related/ukb24983_v2_hg19.BIN_FC10021050.array-combined.PHENO1.glm.logistic.hybrid
323     data/related/ukb24983_v2_hg19.QT_FC10021050.array-combined.PHENO1.glm.linear
323     data/related/ukb24983_v2_hg19.QT_FC21050.array-combined.PHENO1.glm.linear
```

4 of 180 sumstats have 323 lines.

```
$ cat 4_wc_l.tsv | awk '(NR>1 && $1 != 692){print $2}' | sed -e "s/.linear//g" | sed -e "s/.logistic.hybrid//g" | while read f ; do grep "Skipping --glm regression on phenotype 'PHENO1' since # samples <= #"  ${f%.PHENO1.glm}.log ; done | wc -l
4
```

Those were due to

```
Skipping --glm regression on phenotype 'PHENO1' since # samples <= # predictor columns.
```

Then, we identify 227 (= 407 - 180) (pop, GBE_ID) pairs.

```
bash 5_identify.227.sh
```

Check the error messages:

```{bash}
cat 5_identify.227.tsv | awk '(NR>1){print $NF}' | while read f ; do grep "Skipping" $f | uniq  ; done | sort -u
--glm: Skipping case/control phenotype 'PHENO1' since all samples are cases.
--glm: Skipping case/control phenotype 'PHENO1' since all samples are controls.
--glm: Skipping constant quantitative phenotype 'PHENO1'.
Note: Skipping --glm regression on phenotype 'PHENO1' since # samples <= #
```

Check the numbers

```
cat 5_identify.227.tsv | awk '(NR>1){print $NF}' | while read f ; do grep "Skipping case/control phenotype 'PHENO1' since all samples are cases" $f | uniq  ; done | wc -l
2
```

```
cat 5_identify.227.tsv | awk '(NR>1){print $NF}' | while read f ; do grep "Skipping case/control phenotype 'PHENO1' since all samples are controls" $f | uniq  ; done | wc -l
196
```

```
cat 5_identify.227.tsv | awk '(NR>1){print $NF}' | while read f ; do grep "Skipping constant quantitative phenotype 'PHENO1'" $f | uniq  ; done | wc -l
17
```

```
cat 5_identify.227.tsv | awk '(NR>1){print $NF}' | while read f ; do grep "Skipping --glm regression on phenotype 'PHENO1' since # samples <=" $f | uniq ; done | wc -l
12
```

- 407 (pop, GBE_ID)
  - 180 have some results
    - 176 are finished
    - 4 have 322 variants (369 variants are missing)
      - `Skipping --glm regression on phenotype 'PHENO1' since # samples <= # predictor columns.`
  - 227 were skipped
    - 2: all cases
    - 196: all controls
    - 17: constant quantitative phe
    - 12: `# samples <= # predictor columns`

We then combine the results

```{bash}
seq $(cat 6_map.tsv | egrep -v '^#' | wc -l ) > 7_merge.job.lst

sbatch -p mrivas,normal --nodes=1 --mem=4000 --cores=1 --time=1:00:00 --job-name=g691 --output=logs/g691.%A_%a.out --error=logs/g691.%A_%a.err --array=1-180 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 7_merge.sh 7_merge.job.lst 1
```

The summary statistics of the 691 variants were now integrated.
