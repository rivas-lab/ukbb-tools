# Finalizing the GWAS of the 402 variants

## step 1 simple check with `bash` scripts and commands

### step 1-1, copy the files to `/scratch`

- original: `/oak/stanford/groups/mrivas/users/guhan/sandbox/variants_402`
- copy: `/scratch/groups/mrivas/users/ytanigaw/20200703_gwas_merge_bkup/variants_402`

We work on the copy

### step 1-2, check the files in the 402 variants dir

```{bash}
dir402=/scratch/groups/mrivas/users/ytanigaw/20200703_gwas_merge_bkup/variants_402
find ${dir402} -type f -name "*.glm.*" | while read f ; do wc $f ; done | awk '$1 != 403'
```

this file was empty:

- `white_british/ukb24983_v2_hg19.2401-3562.array-combined.HC976.glm.logistic.hybrid`

I deleted this file

### step 1-3, check if we have one file per (GBE_ID, pop)

it is not the case

```
find ${dir402}/white_british -type f -name "*.glm.*" | while read f ; do basename $f ; done | awk -v FS='.' '{print $2}' | sort | uniq -c
    1338 1
    1238 2
     483 3
     398 4
```

### step 1-4, check if we have results for 402 variants

```{bash}
{
    echo "#GBE_ID population f402"
    find ${dir402} -type f -name "*.glm.*" | while read f ; do
        pop=$(basename $(dirname $f))
        GBE_ID=$(basename $f | awk -v FS='.' '{print $4}')
        echo $GBE_ID $pop $f
    done
} | tr " " "\t" > 20200703-gwas-402.1.check.tsv
```

## step 2 identify the list of missing phenotypes for the 402 variant GWAS

- [`20200703-gwas-402.2.check.ipynb`](20200703-gwas-402.2.check.ipynb)
- From `gwas-current-gz-wc.20200629-131527.tsv`, we selected (pop, GBE_ID) where having 402 variants would complete the GWAS analysis.
- We checked if we have such sumstats by looking at [`20200703-gwas-402.1.check.tsv`](20200703-gwas-402.1.check.tsv)
- The results are written to [`20200703-gwas-402.2.still-missing-402.tsv`](20200703-gwas-402.2.still-missing-402.tsv)

## Step 3. run the GWAS for 345 (pop, GBE_ID) pairs

- [`20200703-gwas-402.3.gwas.sh`](20200703-gwas-402.3.gwas.sh)
- We originally used master phe file and tested this script with `e_asian` and `HC411`.
  - Then we found that our master phe file has yet another issue. Some phenotypes are missing.
  - That phenotype is not accessible from `ukbb-query_master_phe_info.sh`
- We decided to use the single-trait phe file and updated the script.

```{bash}
cat 20200703-gwas-402.2.still-missing-402.tsv | egrep -v '^#' | while read GBE_ID population ; do echo $GBE_ID $population; bash 20200703-gwas-402.3.gwas.sh $GBE_ID $population; done
```

### error messages

#### e_asian HC411

```
Error: --pheno-name: 'HC411' does not appear in header line of
/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20200408.phe.
End time: Fri Jul  3 14:32:15 2020
```

## Step 4. check again

```
dir402=/scratch/groups/mrivas/users/ytanigaw/20200703_gwas_merge_bkup/variants_402
find ${dir402} -type f -name "*.glm.*" | sort > 20200703-gwas-402.4.check.lst
```

## Step 5. create a mapping table between the 402 sumstats and 1,080,567 sumstats

## Step 6. merge them

### test



### job

```{bash}
[ytanigaw@sh02-16n04 ~/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas402]$ sbatch -p mrivas,normal --nodes=1 --mem=4000 --cores=1 --time=4:00:00 --job-name=g402 --output=logs/g402.%A_%a.out --error=logs/g402.%A_%a.err --array=1-950 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 6_merge.sh 6_merge.job.lst 15
Submitted batch job 3556936
```
