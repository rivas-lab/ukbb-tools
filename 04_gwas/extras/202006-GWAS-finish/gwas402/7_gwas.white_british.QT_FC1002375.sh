#!/bin/bash
set -beEuo pipefail

ml load plink2/20200314

memory=30000
cpus=6

GBE_ID=QT_FC1002375
pop=white_british

extract_both_arrays=/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/04_gwas/missing_vars.txt
out_dir=data

get_phe_file () {
    local GBE_ID=$1
    cat /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/05_gbe/phenotype_info.tsv \
    | awk -v GBE_ID="${GBE_ID}" '($1 == GBE_ID){print $NF}'
}

get_plink_suffix () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo glm.linear
    else
        echo glm.logistic.hybrid
    fi
}

get_covar_PCs () {
    local pop=$1

        if [ "${pop}" == "others" ] || [ "${pop}" == "related" ] ; then
        echo Global_PC1-Global_PC18
    else
        echo PC1-PC10
    fi
}

if [ ! -d ${out_dir}/${pop} ] ; then mkdir -p ${out_dir}/${pop} ; fi
glm_suffix=$(get_plink_suffix ${GBE_ID})

plink2 \
  --pfile /oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/pgen/ukb24983_cal_hla_cnv vzs \
  --chr 1-22,X,XY,Y,MT \
  --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
  --covar-name age sex Array $(get_covar_PCs ${pop}) N_CNV LEN_CNV \
  --covar-variance-standardize \
  --extract ${extract_both_arrays} \
  --glm skip firth-fallback hide-covar omit-ref no-x-sex \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe \
  --memory ${memory} \
  --out ${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.one_array \
  --pheno $( get_phe_file ${GBE_ID} ) \
  --pheno-quantile-normalize \
  --threads ${cpus} \
  --vif 100000000

s1=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.one_array.PHENO1.${glm_suffix}
s2=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.both_arrays.PHENO1.${glm_suffix}
s0=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.PHENO1.${glm_suffix}

l1=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.one_array.log
l2=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.both_arrays.log
l0=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.log

cat <(echo "#$(basename $l1)") $l1 <(echo "#$(basename $l2)") $l2 > $l0
rm $l1 $l2

if [ -f $s1 ] && [ -f $s2 ] ; then
  cat $s1 <(cat $s2 | egrep -v '^#') > $s0
  rm $s1 $s2
elif [ -f $s1 ] ; then
  mv $s1 $s0
elif [ -f $s2 ] ; then
  mv $s2 $s0
fi

echo $s0
echo $l0
