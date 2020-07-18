# UKB meta-analysis

https://github.com/rivas-lab/ukbb-tools/issues/22


```
cat ../../04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200704-155715.annotated.tsv | egrep -v '^#' | cut -f1 | sort | comm -12 /dev/stdin <(bash 2_list_GBE_IDs.sh | sort) > 2_list_GBE_IDs.$(date +%Y%m%d-%H%M%S).lst
```

`2_list_GBE_IDs.20200705-114626.lst`

```
sbatch -p mrivas,normal,owners --nodes=1 --mem=4000 --cores=1 --time=30:00 --job-name=metal --output=logs/metal.%A_%a.out --error=logs/metal.%A_%a.err --array=1-933 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 1_metal.sh 2_list_GBE_IDs.20200705-114626.lst 2
Submitted batch job 3618473
```

Version 2020/7/17

https://github.com/rivas-lab/ukbb-tools/issues/22
