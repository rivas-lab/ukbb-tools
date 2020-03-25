#!/bin/bash
set -beEuo pipefail

chr=$1
if [ $# -gt 1 ] ; then cpu=$2 ; else cpu=4 ; fi
if [ $# -gt 2 ] ; then mem=$3 ; else mem=40000 ; fi

bgen2pgen_conv () {
    local c=$1
    local cpu=$2
    local mem=$3

    if [ "${c}" == "X" ] ; then
        s_cnt="@@@"
    elif [ "${c}" == "XY" ] ; then
        s_cnt="@@@"
    else
        s_cnt="487297"
    fi

    local plink_common_opts=" --threads ${cpu} --memory ${mem} "
    local bgen_dir="/scratch/groups/mrivas/ukbb24983/hap/download"
    local sample="${bgen_dir}/ukb24983_hap_chr${c}_v2_s${s_cnt}.sample"    
    local scratch_dir="/scratch/groups/mrivas/ukbb24983/hap/pgen"
    local oak_dir="/oak/stanford/groups/mrivas/ukbb/24983/hap/pgen"
    local base="ukb_hap_chr${c}_v2"
    local oak_f="${oak_dir}/${base}"
    local scr_f="${scratch_dir}/${base}"
    local bgen=${bgen_dir}/${base}.bgen
    
    if [ ! -f ${oak_f}.pgen.log ] ; then
        plink2 ${plink_common_opts} --out ${scr_f} \
        --bgen "${bgen}" ref-first \
        --sample "${sample}" \
        --make-pgen vzs \
        --oxford-single-chr ${c}
        mv ${scr_f}.log ${scr_f}.pgen.log
        for ext in pvar.zst psam pgen.log ; do
           cp ${scr_f}.${ext} ${oak_f}.${ext}
        done
        ln -sf ${scr_f}.pgen ${oak_f}.pgen
    fi
    
    if [ ! -f ${oak_f}.bed.log ] ; then
        plink2 ${plink_common_opts} --out ${scr_f} \
        --pgen ${scr_f}.pgen --psam ${oak_f}.psam --pvar ${oak_f}.pvar.zst --make-bed      
        mv ${scr_f}.log ${scr_f}.bed.log
        for ext in bim fam bed.log ; do
           cp ${scr_f}.${ext} ${oak_f}.${ext}
        done
        ln -sf ${scr_f}.bed ${oak_f}.bed
    fi
}

bgen2pgen_conv $chr $cpu $mem
