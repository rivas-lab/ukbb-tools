#!/bin/bash
set -beEuo pipefail

find_col_idx_by_field_id () {
    local tab_file=$1
    local field_id=$2

    local tab_columns_file=${tab_file}.columns

    grep -n ${field_id} ${tab_columns_file} \
    | awk -v FS=':' '{print $1}' | tr "\n" "," | rev | cut -c2- | rev
}

extract_cols () {
    local tab_file=$1
    local cols=$2

    cat ${tab_file} \
    | cut -f "1,${cols}" \
    | sed -e "s/f\\.//g" \
    | sed -e "s/^eid/IID/g"
}

find_tab_file () {
    local basket_id=$1
    local table_id=$2

    local tbl_oak="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/${basket_id}/${table_id}/download/ukb${table_id}.tab"

    echo ${tbl_oak}
}

find_scratch_copy () {
    local oak_file=$1

    local oak="/oak/stanford/groups/mrivas/"
    local scratch="/scratch/groups/mrivas/"

    local scratch_file=$( echo ${oak_file} | sed -e "s%${oak}%${scratch}%g" )


    if [ -f ${scratch_file} ] ; then
        echo ${scratch_file}
    elif [ -f ${oak_file} ] ; then
        echo ${oak_file}
    else
        echo "${oak_file} does not exist!" >&2
        exit 1
    fi
}
