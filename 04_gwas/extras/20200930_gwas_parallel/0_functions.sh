#!/bin/bash
set -beEuo pipefail

source "../functions.sh"


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
