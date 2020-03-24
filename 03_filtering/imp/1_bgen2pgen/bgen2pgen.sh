#!/bin/bash
set -beEuo pipefail

chr=$1
if [ $# -gt 1 ] ; then cpu=$2 ; else cpu=10 ; fi
if [ $# -gt 2 ] ; then mem=$3 ; else mem=150000 ; fi

imp_conv () {
    local c=$1
    local cpu=$2
    local mem=$3

    local plink_common_opts=" --threads ${cpu} --memory ${mem} "
    local bgen_dir="/oak/stanford/projects/ukbb/genotypes/EGAD00010001474"
    local scratch_dir="/scratch/groups/mrivas/ukbb/24983/imp/pgen"
    local oak_dir="/oak/stanford/groups/mrivas/ukbb/24983/imp/pgen"
    local basename="ukb24983_imp_chr${c}_v3"
    local oak_f="${oak_dir}/ukb24983_imp_chr${c}_v3"
    local scr_f="${scratch_dir}/ukb24983_imp_chr${c}_v3"

    if [ "${c}" == "X" ] ; then
        s_cnt="486666"
    elif [ "${c}" == "XY" ] ; then
        s_cnt="486352"
    else
        s_cnt="487317"
    fi

    if [ ! -f ${oak_f}.pgen.log ] ; then
        plink2 ${plink_common_opts} --out ${scr_f} \
        --bgen "${bgen_dir}/ukb_imp_chr${c}_v3.bgen" ref-first \
        --sample "${oak_dir}/ukb24983_imp_chr${c}_v3_s${s_cnt}.sample" \
        --make-pgen vzs
        mv ${scr_f}.log ${oak_f}.pgen.log
        #for ext in pvar.zst psam ; do
        #    mv ${scr_f}.${ext} ${oak_f}.${ext}
        #done
        ln -sf ${scr_f}.pgen ${oak_f}.pgen
    fi
    
    if [ ! -f ${oak_f}.bed.log ] ; then
        plink2 ${plink_common_opts} --out ${scr_f} \
        --pgen ${scr_f}.pgen --psam ${oak_f}.psam --pvar ${oak_f}.pvar.zst --make-bed
        mv ${scr_f}.log ${oak_f}.bed.log
        #for ext in bim fam ; do
        #    mv ${scr_f}.${ext} ${oak_f}.${ext}
        #done
        ln -sf ${scr_f}.bed ${oak_f}.pgen
    fi
}

imp_conv $chr $cpu $mem

