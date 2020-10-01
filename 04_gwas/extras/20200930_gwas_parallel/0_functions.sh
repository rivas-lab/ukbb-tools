#!/bin/bash
set -beEuo pipefail

load_plink2 () {
    # We installed PLINK2 software as a software module in our HPC system.
    # We call `ml load plink2` for the specified version. By doing this, 
    # it updates the PATHs so that we can execute plink2 software.

    plink2_version=$1

    check_avx2_flag=$(cat /proc/cpuinfo | grep flags | grep -i avx2 | uniq | cat /dev/stdin <(echo avx2) | grep -i avx2 | wc -l)

    if [ ${check_avx2_flag} -eq 2 ] ; then
        ml load plink2/${plink2_version}
    else
        ml load plink2/${plink2_version}-non-AVX2
    fi
}

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
    cat /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/05_gbe/phenotype_info.tsv \
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

get_covar_names () {
    local pop=$1
    local batch_idx=$2
    local n_batch_one_array=$3

    if [ ${batch_idx} -le ${n_batch_one_array} ] ; then
        echo "age,sex,$(get_covar_PCs ${pop}),N_CNV,LEN_CNV"
    else
        echo "age,sex,Array,$(get_covar_PCs ${pop}),N_CNV,LEN_CNV"
    fi
}
