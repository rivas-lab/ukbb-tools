#!/bin/bash
set -beEuo pipefail

ml load plink2/20200314

memory=30000
cpus=6

# GBE_ID=$1
# pop=$2

GBE_ID=INI22149
pop=others
extract=missing.369.lst
out_dir=data

get_phe_file () {
    local GBE_ID=$1
    cat /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/05_gbe/phenotype_info.tsv \
    | awk -v GBE_ID="${GBE_ID}" '($1 == GBE_ID){print $NF}'
}

if [ ! -d ${out_dir}/${pop} ] ; then mkdir -p ${out_dir}/${pop} ; fi

plink2 \
  --pfile /oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/pgen/ukb24983_cal_hla_cnv vzs \
  --chr 1-22,X,XY,Y,MT \
  --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
  --covar-name age sex Array Global_PC1-Global_PC18 N_CNV LEN_CNV \
  --covar-variance-standardize \
  --extract ${extract} \
  --glm skip firth-fallback hide-covar omit-ref no-x-sex \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe \
  --memory ${memory} \
  --out ${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined \
  --pheno $( get_phe_file ${GBE_ID} ) \
  --pheno-quantile-normalize \
  --threads ${cpus} \
  --vif 100000000

exit 0
