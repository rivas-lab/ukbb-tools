#!/bin/bash
set -beEuo pipefail

[[ ${GRM_CAL_BFILE:-}  -eq 1 ]]  && return || readonly GRM_CAL_BFILE="/scratch/groups/mrivas/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2_hg19"
#[[ ${GRM_CAL_BFILE:-}  -eq 1 ]]  && return || readonly GRM_CAL_BFILE="/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19"
[[ ${GRM_CAL_DIR:-}  -eq 1 ]]    && return || readonly GRM_CAL_DIR="/oak/stanford/groups/mrivas/ukbb24983/cal/grm"

get_keep () {
    local pop=$1
    echo "/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_${pop}.phe"
}

gcta_grm_part () {    
    local threads=$1
    local pop=$2
    local batch_num=$3
    local n_parts=$4 # 100

    local bfile=${GRM_CAL_BFILE}
    local keep=$(get_keep $pop)
    local out="${GRM_CAL_DIR}/part/ukb24983_cal_cALL_v2_hg19.${pop}"

    gcta64 \
        --make-grm-part $n_parts $batch_num \
        --bfile ${bfile} \
        --keep ${keep} \
        --maf 0.01 \
        --thread-num ${threads} \
        --out ${out}
}

gcta_grm_combine () {
    local pop=$1
    local n_parts=$2 # 100
    local part="${GRM_CAL_DIR}/part/ukb24983_cal_cALL_v2_hg19.${pop}"
    local out="${GRM_CAL_DIR}/ukb24983_cal_cALL_v2_hg19.${pop}"

    for ext in log grm.id grm.bin grm.N.bin ; do
        for i in $(seq -w $n_parts) ; do
            file_part="${part}.part_${n_parts}_${i}.${ext}"
            if [ -f ${file_part} ] ; then
                cat ${file_part}
            else
                echo [warning] missing_file: ${file_part} >&2
            fi
        done > ${out}.${ext}
    done
}

# gcta_grm_part_find_unfinished () {
#     comm --nocheck-order -3 <(seq 100) <(cat logs/GRM.*_*.err | grep array-end | awk -v FS='=' '{print $NF}'| sort -nu )
# }

# gcta_grm_part_find_intermediate_files () {
#     local res_dir="/oak/stanford/groups/mrivas/ukbb24983/cal/grm/part/"
#     find_unfinished | while read i ; do 
#         find ${res_dir} -name "ukb24983_cal_cALL_v2_hg19.white_british.part_100_$(printf "%03d" $i).*" 
#     done
# }
