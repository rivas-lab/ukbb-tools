#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}

if [ $# -gt 0 ] ; then batch_idx=$1 ; fi
# if [ $# -gt 1 ] ; then GBE_ID=$2 ; else GBE_ID=HC382 ; fi
# if [ $# -gt 2 ] ; then plink2_version=$3 ; else plink2_version=20200725 ; fi

memory=7000
cpus=2

GBE_ID=HC382
plink2_version=20200727
pop=white_british
n_batch=100
# pfile="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv"
pfile="/scratch/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv"
one_array="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/one_array_variants.txt"
both_arrays="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/both_arrays_variants.txt"
out_dir="/oak/stanford/groups/mrivas/users/ytanigaw/sandbox/20200725-plink2-firth"

#########################

compute_n_batch_one_array () {
    n_batch=$1
    pfile=$2
    one_array=$3

    # check if this batch idx correspond to the "one array" list or "both arrays" list
    one_array_n=$(cat ${one_array} | wc -l)
    pvar_n_vars=$(zstdcat ${pfile}.pvar.zst | egrep -v '^#' | wc -l)
    perl -e "print(  int(  0.5 + ${n_batch} * ${one_array_n} / ${pvar_n_vars}  )  )"
}

show_var_list () {
    batch_idx=$1
    n_batch=$2
    n_batch_one_array=$3
    one_array=$4
    both_arrays=$5

    if [ ${batch_idx} -le ${n_batch_one_array} ] ; then
        var_list=${one_array}
        n_batch_this=${n_batch_one_array}
        batch_idx_offset=0
    else
        var_list=${both_arrays}
        n_batch_this=$(perl -e "print(  ${n_batch} - ${n_batch_one_array}  )")
        batch_idx_offset=${n_batch_one_array}
    fi

    var_n=$(cat ${var_list} | egrep -v '^#' | wc -l)
    batch_size=$(perl -e "print(  int((${var_n}-1)/${n_batch_this}) + 1  )")
    idx_s=$(perl -e "print(  ${batch_size} *  (${batch_idx} - ${batch_idx_offset} - 1) + 1  )")
    idx_e=$(perl -e "print(  ${batch_size} *  (${batch_idx} - ${batch_idx_offset}) )")

    cat ${var_list} | egrep -v '^#' | awk -v idx_s=${idx_s} -v idx_e=${idx_e} '(idx_s <= NR && NR <= idx_e)'
}

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

get_covar_name () {
    local pop=$1
    local batch_idx=$2
    local n_batch_one_array=$3

    if [ ${batch_idx} -le ${n_batch_one_array} ] ; then
        echo "age sex $(get_covar_PCs ${pop}) N_CNV LEN_CNV"
    else
        echo "age sex Array $(get_covar_PCs ${pop}) N_CNV LEN_CNV"
    fi
}


#########################

check_avx2_flag=$(cat /proc/cpuinfo | grep flags | grep -i avx2 | uniq | cat /dev/stdin <(echo avx2) | grep -i avx2 | wc -l)

if [ ${check_avx2_flag} -eq 2 ] ; then
    ml load plink2/${plink2_version}
else
    ml load plink2/${plink2_version}-non-AVX2
fi

# We installed PLINK2 software as a software module in our HPC system.
# This `ml load plink2/20200409` updates the PATHs so that we can execute plink2 software.

# if [ ! -d ${out_dir}/${pop} ] ; then mkdir -p ${out_dir}/${pop} ; fi

if [ ! -d ${out_dir}/cc-residualize ] ; then mkdir -p ${out_dir}/cc-residualize ; fi
if [ ! -d ${out_dir}/firth-residualize ] ; then mkdir -p ${out_dir}/firth-residualize ; fi
if [ ! -d ${out_dir}/firth-default ] ;     then mkdir -p ${out_dir}/firth-default ; fi

n_batch_one_array=$(compute_n_batch_one_array ${n_batch} ${pfile} ${one_array})
glm_suffix=$(get_plink_suffix ${GBE_ID})

# cc-residualize

show_var_list ${batch_idx} ${n_batch} ${n_batch_one_array} ${one_array} ${both_arrays} \
| plink2 \
  --memory ${memory} \
  --threads ${cpus} \
  --pfile ${pfile} vzs \
  --chr 1-22,X,XY,Y,MT \
  --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
  --covar-name $( get_covar_name ${pop} ${batch_idx} ${n_batch_one_array} ) \
  --extract /dev/stdin \
  --glm skip-invalid-pheno firth-fallback firth-residualize cc-residualize hide-covar omit-ref no-x-sex \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe \
  --out ${out_dir}/cc-residualize/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx} \
  --pheno $( get_phe_file ${GBE_ID} ) \
  --covar-variance-standardize \
  --pheno-quantile-normalize \
  --vif 100000000

mv \
${out_dir}/cc-residualize/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.PHENO1.${glm_suffix} \
${out_dir}/cc-residualize/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.${glm_suffix}

# firth-residualize

show_var_list ${batch_idx} ${n_batch} ${n_batch_one_array} ${one_array} ${both_arrays} \
| plink2 \
  --memory ${memory} \
  --threads ${cpus} \
  --pfile ${pfile} vzs \
  --chr 1-22,X,XY,Y,MT \
  --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
  --covar-name $( get_covar_name ${pop} ${batch_idx} ${n_batch_one_array} ) \
  --extract /dev/stdin \
  --glm skip-invalid-pheno firth-fallback firth-residualize hide-covar omit-ref no-x-sex \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe \
  --out ${out_dir}/firth-residualize/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx} \
  --pheno $( get_phe_file ${GBE_ID} ) \
  --covar-variance-standardize \
  --pheno-quantile-normalize \
  --vif 100000000

mv \
${out_dir}/firth-residualize/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.PHENO1.${glm_suffix} \
${out_dir}/firth-residualize/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.${glm_suffix}

# firth default

show_var_list ${batch_idx} ${n_batch} ${n_batch_one_array} ${one_array} ${both_arrays} \
| plink2 \
  --memory ${memory} \
  --threads ${cpus} \
  --pfile ${pfile} vzs \
  --chr 1-22,X,XY,Y,MT \
  --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
  --covar-name $( get_covar_name ${pop} ${batch_idx} ${n_batch_one_array} ) \
  --extract /dev/stdin \
  --glm skip-invalid-pheno firth-fallback hide-covar omit-ref no-x-sex \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe \
  --out ${out_dir}/firth-default/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx} \
  --pheno $( get_phe_file ${GBE_ID} ) \
  --covar-variance-standardize \
  --pheno-quantile-normalize \
  --vif 100000000

mv \
${out_dir}/firth-default/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.PHENO1.${glm_suffix} \
${out_dir}/firth-default/ukb24983_v2_hg19.${GBE_ID}.batch${batch_idx}.${glm_suffix}
