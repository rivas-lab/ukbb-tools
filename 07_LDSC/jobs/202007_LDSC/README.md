# LDSC 2020/07

## apply LDSC munge

```{bash}
bash 1_generate_input_list.sh | tee 1_LDSC_munge.$(date +%Y%m%d-%H%M%S).job.lst | tee /dev/stderr | wc -l

ml load resbatch

sbatch -p mrivas,normal,owners --time=3:00:00 --mem=8000 --nodes=1 --cores=1 --job-name=munge --output=logs/munge.%A_%a.out --error=logs/munge.%A_%a.err --array=1-905 $parallel_sbatch_sh 1_LDSC_munge.sh 1_LDSC_munge.20200717-210250.job.lst 3
```

- https://github.com/rivas-lab/ukbb-tools/issues/26
### Check the completion of the scripts

```
find /oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc -type f -name "*.gz" | wc
```

```
find logs/ -name "*.err" | while read f ; do cat $f | grep array-end | awk -v FS='=' '{print $NF}'; done | sort | comm -3 <(seq 955 | sort) /dev/stdin | sort -n | tr '\n' ','
```

## apply LDSC h2

```
pop="white_british"
cat ../../../04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200704-155715.combined.tsv | awk -v FS='\t' -v wc_l=1080969 -v pop=${pop} '($5 == wc_l && $2 == pop){print $1}' > 2_ldsc_h2.${pop}.$(date +%Y%m%d-%H%M%S).job.lst
```

```{bash}
sbatch -p mrivas,normal,owners --nodes=1 --mem=8000 --cores=1 --time=1:00:00 --job-name=UKBh2 --output=logs/UKBh2.%A_%a.out --error=logs/UKBh2.%A_%a.err --array=1-897 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 2_ldsc_h2.sh 2_ldsc_h2.white_british.20200705-212236.job.lst 4
Submitted batch job 3656875
```
