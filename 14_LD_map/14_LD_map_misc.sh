#!/bin/bash
set -beEuo pipefail

compute_ld_bool_wrapper () {
    # run plink2 --indep-pairwise
    local bfile=$1
    local out_prefix=$2
    local cores=$3
    local memory=$4
    if [ $# -gt 4 ] ; then local keep_file=$5  ; else local keep_file="" ; fi
    if [ $# -gt 5 ] ; then local plink_opts=$6 ; else local plink_opts="" ; fi
    local plink_opts_full="$(echo ${plink_opts} | tr ":" " ") --memory ${memory} --threads ${cores}"

    if [ ! -d $(dirname ${out_prefix}) ] ; then mkdir -p $(dirname ${out_prefix}) ; fi
    
    if [ ! -f ${out_prefix}.bool.prune.in ] && [ ! -f ${out_prefix}.bool.prune.out ] ; then
        plink2 ${plink_opts_full} \
            --pfile ${bfile} vzs $([ "${keep_file}" == "" ] && echo "" || echo "--keep ${keep_file}") \
            --allow-extra-chr \
            --indep-pairwise 50 5 0.5 \
            --out ${out_prefix}.bool
    fi   
}

compute_ld_r2_wrapper () {
    # run plink1.9 --r2
    local bfile=$1
    local out_prefix=$2
    local cores=$3
    local memory=$4
    if [ $# -gt 4 ] ; then local keep_file=$5  ; else local keep_file="" ; fi
    if [ $# -gt 5 ] ; then local plink_opts=$6 ; else local plink_opts="" ; fi
    local plink_opts_full="$(echo ${plink_opts} | tr ":" " ") --memory ${memory} --threads ${cores}"

    if [ ! -d $(dirname ${out_prefix}) ] ; then mkdir -p $(dirname ${out_prefix}) ; fi
    
    if [ ! -f ${out_prefix}.ld_map.ld.gz ] ; then
        plink ${plink_opts_full} \
            --bfile ${bfile} $([ "${keep_file}" == "" ] && echo "" || echo "--keep ${keep_file}") \
            --allow-extra-chr \
            --ld-window-kb 1000 --ld-window-r2 0.1 --r2 gz \
            --ld-window 100000000 \
            --out ${out_prefix}.ld_map    
    fi

    if [ ! -f ${out_prefix}.ld_map.tsv.gz ] ; then
        zcat ${out_prefix}.ld_map.ld.gz \
            | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7}' \
            | sed -e "s/^CHR_A/#CHR_A/g" \
            | bgzip -l9 -@${cores} > ${out_prefix}.ld_map.tsv.gz
    fi

    if [ ! -f ${out_prefix}.ld_map.tsv.gz.tbi ] ; then
        tabix -c '#' -s 1 -b 2 -e 5 ${out_prefix}.ld_map.tsv.gz 
    fi
}

compute_ld_map_wrapper () {
    local bfile=$1
    local out_prefix=$2
    local cores=$3
    local memory=$4
    if [ $# -gt 4 ] ; then local keep_file=$5  ; else local keep_file="" ; fi
    if [ $# -gt 5 ] ; then local plink_opts=$6 ; else local plink_opts="" ; fi
    local plink_opts_full="$(echo ${plink_opts} | tr ":" " ") --memory ${memory} --threads ${cores}"

    if [ ! -d $(dirname ${out_prefix}) ] ; then mkdir -p $(dirname ${out_prefix}) ; fi

    compute_ld_bool_wrapper $@

    compute_ld_r2_wrapper $@
}

check_ldmap_out () {
    local basename=$1

    for ext in "bool.prune.in" "bool.prune.out" "ld_map.ld.gz" "ld_map.tsv.gz" "ld_map.tsv.gz.tbi" ; do
        local file="${basename}.${ext}"
        if [ ! -f ${file} ] ; then
            echo "$(basename ${basename}) ${file}" | tr " " "\t"
        fi
    done
}
