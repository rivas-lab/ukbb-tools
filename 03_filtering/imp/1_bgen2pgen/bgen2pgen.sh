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
    local scratch_dir="/scratch/groups/mrivas/ukbb/24983/imp/pgen"
    local out_prefix="ukb24983_imp_chr${c}_v3"

    if [ ! -f ${out_prefix}.pgen.log ] && [ ! -f ${out_prefix}-temporary.psam ] ; then
        plink2 ${plink_common_opts} --out ${scratch_dir}/${out_prefix} \
        --bgen $(readlink -f ukb_imp_chr${c}_v3.bgen) ref-first \
        --sample $(readlink -f $(find $(pwd) -name "ukb24983_imp_chr${c}_v3_s*.sample")) \
        --make-pgen vzs
        mv ${scratch_dir}/${out_prefix}.log ${out_prefix}.pgen.log
        for ext in pvar.zst psam ; do
            mv ${scratch_dir}/${out_prefix}.${ext} ${out_prefix}.${ext}
        done
        ln -s ${scratch_dir}/${out_prefix}.pgen .
    fi
    
    if [ ! -f ${out_prefix}.bed.log ]  && [ ! -f ${out_prefix}-temporary.bim ]; then
        plink2 ${plink_common_opts} --out ${scratch_dir}/${out_prefix} \
        --pfile ${out_prefix} vzs --make-bed
        mv ${scratch_dir}/${out_prefix}.log ${out_prefix}.bed.log
        for ext in bim fam ; do
            mv ${scratch_dir}/${out_prefix}.${ext} ${out_prefix}.${ext}
        done
        ln -s ${scratch_dir}/${out_prefix}.bed .
    fi
}

imp_conv $chr $cpu $mem

