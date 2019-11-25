#!/bin/bash
set -beEuo pipefail

gcta_grm_part () {
    local batch_num=$1
    local threads=$2
    local pop="white_british"

    local bfile="/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19"
    local keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_${pop}.phe"
    local out="/oak/stanford/groups/mrivas/ukbb24983/cal/grm/part/ukb24983_cal_cALL_v2_hg19.${pop}"

    gcta64 \
        --make-grm-part 1000 $batch_num \
        --bfile ${bfile} \
        --keep ${keep} \
        --maf 0.01 \
        --thread-num ${threads} \
        --out ${out}
}

gcta_grm_combine () {
    local pop=$1
    local part="/oak/stanford/groups/mrivas/ukbb24983/cal/grm/part/ukb24983_cal_cALL_v2_hg19.${pop}"
    local out="/oak/stanford/groups/mrivas/ukbb24983/cal/grm/ukb24983_cal_cALL_v2_hg19.${pop}"

    local n_parts=1000

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

