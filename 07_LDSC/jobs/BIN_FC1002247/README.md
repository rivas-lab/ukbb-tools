```
find /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/phe/HC0.phe /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/extras/highconfidenceqc/white_british/ -maxdepth 1 -type f -name "*.gz" > HC.sumstats.lst

cat HC.sumstats.lst | parallel -k 'echo {} /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/9796/24611/white_british/ukb24983_v2_hg19.BIN_FC1002247.genotyped.PHENO1.glm.logistic.hybrid.gz' | tr " " "\t" > jobs.tsv

echo "#p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se" | tr " " "\t" ; find /oak/stanford/groups/mrivas/dev-ukbb-tools/ldsc/rg/BIN_FC1002247/ -maxdepth 1 -type f -name "*.log" | parallel 'cat {} | tail -n4 | head -n1 ' | sed -e "s%/oak/stanford/groups/mrivas/dev-ukbb-tools/ldsc/munged/BIN_FC1002247/ukb24983_v2_hg19.%%g" | sed -e "s%.genotyped.sumstats.gz%%g" | sort -k1V | awk -v OFS='\t' '{ for (i = 1; i <= 12; i++) { printf $i "\t" } print "\n"}' | sed -e "s/\t$//g"
```

