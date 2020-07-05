# LDSC 2020/07

## apply LDSC munge

```
cat ~/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200704-155715.annotated.tsv | awk 'NR==1 || $NF == 1080969' > 1_LDSC_munge.$(date +%Y%m%d-%H%M%S).tsv
```

`1_LDSC_munge.20200704-171713.tsv`
20940 files

21 files each, 998 tasks

```
cat 1_LDSC_munge.20200704-171713.tsv | egrep -v '#' | cut -f3  > 1_LDSC_munge.20200704-171713.job.lst

sbatch -p mrivas,normal,owners --nodes=1 --mem=8000 --cores=1 --time=6:00:00 --job-name=munge --output=logs/munge.%A_%a.out --error=logs/munge.%A_%a.err --array=1-998 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 1_LDSC_munge.sh 1_LDSC_munge.20200704-171713.job.lst 21
Submitted batch job 3581032
```