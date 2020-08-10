# UKB meta-analysis

We perform UKB-wide meta-analysis across 7 populations: white British, non-British white, African, South Asian, East Asian, semi-related individuals, and "others" group.

https://github.com/rivas-lab/ukbb-tools/issues/22

```
cat ../../04_gwas/extras/202006-GWAS-finish/gwas-current-gz-wc.20200704-155715.annotated.tsv | egrep -v '^#' | cut -f1 | sort | comm -12 /dev/stdin <(bash 2_list_GBE_IDs.sh | sort) > 2_list_GBE_IDs.$(date +%Y%m%d-%H%M%S).lst
```

`2_list_GBE_IDs.20200705-114626.lst`

```
sbatch -p mrivas,normal,owners --nodes=1 --mem=4000 --cores=1 --time=30:00 --job-name=metal --output=logs/metal.%A_%a.out --error=logs/metal.%A_%a.err --array=1-933 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 1_metal.sh 2_list_GBE_IDs.20200705-114626.lst 2
Submitted batch job 3618473
```

## Version 2020/7/17

For 3,794 phenotypes, we have the results for UKB-wide meta-analysis. This analysis was performed by Yosuke Tanigawa.

https://github.com/rivas-lab/ukbb-tools/issues/22

### `icdinfo` file

We checked the log-files from Metal and identify the list of population actually used in the meta-analysis for each phenotype ([`4_list_metal_inputs.sh`](4_list_metal_inputs.sh)). We looked at the phenotype info file and computed the sum of Ns across the identified populations to tabulate the meta-analysis Ns ([`4_icdinfo.ipynb`](4_icdinfo.ipynb)).
