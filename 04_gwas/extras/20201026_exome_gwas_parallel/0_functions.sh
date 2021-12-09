#!/bin/bash
set -beEuo pipefail

source "../functions.sh"

show_var_list () {
    batch_idx=$1
    n_batch=$2
    pvar=$3

    var_n=$(cat_or_zcat ${pvar} | egrep -v '^#' | wc -l)
    batch_size=$(perl -e "print(  int((${var_n}-1)/${n_batch}) + 1  )")
    idx_s=$(perl -e "print(  ${batch_size} * (${batch_idx} - 1) + 1  )")
    idx_e=$(perl -e "print(  ${batch_size} *  ${batch_idx})")

    cat_or_zcat ${pvar} | egrep -v '^#' | cut -f3 | awk -v idx_s=${idx_s} -v idx_e=${idx_e} '(idx_s <= NR && NR <= idx_e)'
}

get_covar_names () {
    local pop=$1
    echo "age,sex,$(get_covar_PCs ${pop})"
}

######################################
# look-up file
######################################

master_gwas_dump_file () {
    local f=$1
    GBE_ID=$(echo $f | awk -v FS='.' '{print $2}')
    pop=$(basename $(dirname $f))
    cat_or_zcat $f | egrep -v '^#' \
    | awk -v pop=${pop} -v phe=${GBE_ID} -v OFS='\t' '{print pop, phe, $0}'
}

