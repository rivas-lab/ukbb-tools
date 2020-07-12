#!/bin/bash
set -beEuo pipefail

ml load plink2/20200314

memory=30000
cpus=8

batch_idx=${SLURM_ARRAY_TASK_ID:=204}
if [ $# -gt 0 ] ; then batch_idx=$1 ; fi

jobs="1_gwas.jobs.20200709-081630.tsv"

GBE_ID=$( cat ${jobs} | egrep -v '^#' | awk -v nr=${batch_idx} '(NR==nr){print $1}' )
pop=$(    cat ${jobs} | egrep -v '^#' | awk -v nr=${batch_idx} '(NR==nr){print $2}' )

extract_one_array="/oak/stanford/groups/mrivas/ukbb24983/sqc/one_array_variants.txt"
out_dir="data"

##########################
# load the modules
##########################

is_scg=$(echo scg $(hostname) | tr " " "\n" | grep scg | wc -l)
avx2_flag=$( cat /proc/cpuinfo  | grep flags | uniq  | awk -v FS=':' '{print $2}' | tr " " "\n" | cat /dev/stdin <(echo avx2) | grep -i avx2 | wc -l)
if [ "${is_scg}" -eq 2 ] ; then
    export MODULEPATH="$HOME/.modules:/labs/mrivas/.modules:$MODULEPATH"
    ml load htslib python/3.8.2
else
    ml load htslib python/3.6.1
fi
if [ "${avx2_flag}" -eq 2 ] ; then
    ml load plink2/20200314
else
    ml load plink2/20200314-non-AVX2
fi

##########################
# functions
##########################

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

plink2_glm () {
  plink2 \
    --chr 1-22,X,XY,Y,MT \
    --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
    --covar-variance-standardize \
    --glm skip firth-fallback hide-covar omit-ref no-x-sex \
    --pheno-quantile-normalize \
    --vif 100000000 \
    $@
}

##########################
# run plink2
##########################

if [ ! -d ${out_dir}/${pop} ] ; then mkdir -p ${out_dir}/${pop} ; fi

plink2_glm \
  --threads ${cpus} \
  --memory ${memory} \
  --pfile /oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/pgen/ukb24983_cal_hla_cnv vzs \
  --pheno $( get_phe_file ${GBE_ID} ) \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe \
  --covar-name age sex Array $(get_covar_PCs ${pop}) N_CNV LEN_CNV \
  --exclude ${extract_one_array} \
  --out ${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.both_arrays

plink2_glm \
  --threads ${cpus} \
  --memory ${memory} \
  --pfile /oak/stanford/groups/mrivas/private_data/ukbb/24983/array-combined/pgen/ukb24983_cal_hla_cnv vzs \
  --pheno $( get_phe_file ${GBE_ID} ) \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe \
  --covar-name age sex $(get_covar_PCs ${pop}) N_CNV LEN_CNV \
  --extract ${extract_one_array} \
  --out ${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.one_array

##########################
# combine the log files
##########################

l1=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.one_array.log
l2=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.both_arrays.log
l0=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.log

cat <(echo "#$(basename $l1)") $l1 <(echo "#$(basename $l2)") $l2 > $l0
rm $l1 $l2

##########################
# combine the summary statistics
##########################

glm_suffix=$(get_plink_suffix $GBE_ID)

s1=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.one_array.PHENO1.${glm_suffix}
s2=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.both_arrays.PHENO1.${glm_suffix}
s0=${out_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.${glm_suffix}

if [ -f $s1 ] && [ -f $s2 ] ; then
  cat <(cat $s1 | egrep -v '^#') <(cat $s2 | egrep -v '^#') \
  | sort --parallel ${cpus} -k1,1V -k2,2n -k3,3 \
  | cat <(cat $s1 | egrep '^#') /dev/stdin \
  | bgzip -l9 -@${cpus} > ${s0}.gz
  rm $s1 $s2
elif [ -f $s1 ] ; then
  mv $s1 $s0
  bgzip -l9 -@${cpus} ${s0}
elif [ -f $s2 ] ; then
  mv $s2 $s0
  bgzip -l9 -@${cpus} ${s0}
fi

# show the results files
if [ -f ${s0}.gz ] ; then
  readlink -f ${s0}.gz
fi

# show the log files
readlink -f $l0
