#!/bin/bash
set -beEuo pipefail


# input
metal_list_f="metal2plink.lst"

# output
metal_input_pops="metal_input_pops.tsv"

list_input_files () {
    local info_txt=$1
    cat ${info_txt} | egrep -B10 '^# == METAL info file ==' | egrep -v '^#'
}

list_pops () {
    local info_txt=$1
    list_input_files ${info_txt} | rev | awk -v FS='/' -v OFS='\t' '{print $2, $1}' | rev \
    | sed -e 's/ukb24983_v2_hg19.//g' \
    | sed -e 's/.array-combined.glm.logistic.hybrid.gz//g' \
    | sed -e 's/.array-combined.glm.linear.gz//g'
}

export -f list_input_files
export -f list_pops


{
    echo "#GBE_ID pop" | tr ' ' '\t'

    ! cat ${metal_list_f} | sed -e 's/.metal.tsv.gz$/.metal.info.txt/g' | parallel --eta -j+0 -k list_pops
} > ${metal_input_pops}

echo ${metal_input_pops}
