# LDSC 2020/07

## apply LDSC munge

```{bash}
bash 1_generate_input_list.sh | tee 1_LDSC_munge.$(date +%Y%m%d-%H%M%S).job.lst | tee /dev/stderr | wc -l

sbatch -p mrivas,normal,owners --nodes=1 --mem=8000 --cores=1 --time=6:00:00 --job-name=munge --output=logs/munge.%A_%a.out --error=logs/munge.%A_%a.err --array=1-995 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 1_LDSC_munge.sh 1_LDSC_munge.20200705-135957.job.lst 6
```

- 5727 files, 955 * 6

### Check the completion of the scripts

```
find /oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc -type f -name "*.gz" | wc
```

```
find logs/ -name "*.err" | while read f ; do cat $f | grep array-end | awk -v FS='=' '{print $NF}'; done | sort | comm -3 <(seq 955 | sort) /dev/stdin | sort -n | tr '\n' ','
```
