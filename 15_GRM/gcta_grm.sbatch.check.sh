#!/bin/bash
set -beEuo pipefail

find_unfinished () {
    comm --nocheck-order -3 <(seq 1000) <(cat logs/GRM.*_*.err | grep array-end | awk -v FS='=' '{print $NF}'| sort -nu )
}

find_intermediate_files () {
    local res_dir="/oak/stanford/groups/mrivas/ukbb24983/cal/grm/part/"
    find_unfinished | while read i ; do 
        find ${res_dir} -name "ukb24983_cal_cALL_v2_hg19.white_british.part_1000_$(printf "%04d" $i).*" 
    done
}

find_unfinished
#find_intermediate_files

