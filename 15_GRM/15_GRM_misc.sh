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

