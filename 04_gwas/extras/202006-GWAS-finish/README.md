# GWAS freeze, July 2020

This directory contains scripts to finish up the [GWAS freez, July 2020](https://github.com/rivas-lab/ukbb-tools/milestone/1).

## missing files

```{bash}
bash check.missing_pop_GBE.sh
```

This script checks the phenotype info file and see if we have sumstats.

## `wc -l` analysis

### compute `wc -l` for each file

[`gwas-current-gz-wc.sh`](gwas-current-gz-wc.sh) has more details

```
find /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current -name "*.gz" -type l | sort > /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas-current-gz-list.$(date +%Y%m%d-%H%M%S).txt
```

this generated a list of 22139 files in `gwas-current-gz-list.20200704-150033.txt`

--> 23 each in 963 tasks

```
sbatch -p mrivas,normal --nodes=1 --mem=4000 --cores=1 --time=30:00 --job-name=wc --output=logs/wc.%A_%a.out --error=logs/wc.%A_%a.err --array=1-963 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh gwas-current-gz-wc.sh gwas-current-gz-list.20200704-150033.txt 23

Submitted batch job 3577852
```

### post-processing

```
bash gwas-current-gz-wc-cat.sh
rm gwas-current-gz-list.20200704-150033.txt
find gwas-current-gz-wc -type f | while read f ; do rm $f ; done
```

## combine the results

Using [`gwas-wc-post-processing-20200704.ipynb`](gwas-wc-post-processing-20200704.ipynb), we combined the results from the two types of check shown above.

[`gwas-current-gz-wc.20200704-155715.combined.tsv`](gwas-current-gz-wc.20200704-155715.combined.tsv)

```{bash}
$ cat gwas-current-gz-wc.20200704-155715.combined.tsv | awk '$5 < 1059397' | wc
    302    2114   43103

$ cat gwas-current-gz-wc.20200704-155715.combined.tsv | awk '$5 < 1059397' | cut -f2 | sort | uniq -c
     20 african
     96 e_asian
     41 non_british_white
     12 others
     11 related
     14 s_asian
    108 white_british
```

