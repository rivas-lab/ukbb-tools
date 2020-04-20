# LDSC for the GWAS refresh in 2020/04

The example shown here is for the [COVID-19](https://github.com/rivas-lab/covid19) project.
The job submission script can be found in that directory.

## Step 1. LDSC munge wrapper

```{bash}
bash  ~/repos/rivas-lab/ukbb-tools/07_LDSC/jobs/202004_LDSC/1_LDSC_munge.sh INI30150
```

## Step 2. LDSC rg

```
bash  ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ukb_ldsc_rg_helper.sh \
/scratch/groups/mrivas/projects/covid19/GWAS/LDSC/LDSC.covid19_test_result.INI30190.tsv \
/scratch/groups/mrivas/projects/covid19/GWAS/LDSC/ukb24983_v2_hg19.white_british.array-combined.covid19_test_result.glm.logistic.hybrid.sumstats.gz \
/scratch/groups/mrivas/ukbb24983/array_combined/ldsc/current/white_british/ukb24983_v2_hg19.INI30190.array-combined.glm.sumstats.gz
```

