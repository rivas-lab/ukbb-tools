#!/bin/bash
set -beEuo pipefail

ml load plink2 R/3.6 gcc

# Csq="ptv"
if [ $# -gt 0 ] ; then Csq=$1 ; else Csq="" ; fi

pfile="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv"
QCed_vars_f="/oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/ld_indep_20201015/ukb24983_cal_hla_cnv.var_QCed.tsv.gz"
out_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/ld_indep_20201015/plink_output"
pops=('white_british')

plink_wrapper () {
    plink2 --silent --threads 6 --memory 60000 --pfile ${pfile} vzs $@
}

plink_indep_pairwise () {
#     plink_wrapper --allow-extra-chr --indep-pairwise 50 5 0.5 $@
    plink_wrapper --allow-extra-chr --indep-pairwise 1000kb 1 0.5 $@
}

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] || [ "${file%.bgz}.bgz" == "${file}" ] ; then 
        zcat ${file} 
    elif [ "${file%.zst}.zst" == "${file}" ] ; then 
        zstdcat ${file}
    else
        cat ${file}
    fi
}

show_var_lst () {
    Csq=$1
    
    if   [ ${Csq} == "ptv" ] ; then
        cat_or_zcat ${QCed_vars_f} | awk -v FS='\t' -v Csq=${Csq} '($NF == Csq){print $3}'
    elif [ ${Csq} == "pav" ] ; then
        Rscript 6b_LD_indep_var_filter.R ${Csq} 'ptv' /dev/stdout
    elif [ ${Csq} == "pcv" ] ; then
        Rscript 6b_LD_indep_var_filter.R ${Csq} 'ptv,pav' /dev/stdout
    elif [ ${Csq} == "utr" ] ; then
        Rscript 6b_LD_indep_var_filter.R ${Csq} 'ptv,pav,pcv' /dev/stdout
    elif [ ${Csq} == "intron" ] ; then
        Rscript 6b_LD_indep_var_filter.R ${Csq} 'ptv,pav,pcv,utr' /dev/stdout
    elif [ ${Csq} == "others" ] ; then
        Rscript 6b_LD_indep_var_filter.R ${Csq} 'ptv,pav,pcv,utr,intron' /dev/stdout
    else
        echo "unsupported variant consequence type: ${Csq}" >&2
        exit 1
    fi
}

if [ ! -d ${out_d} ] ; then mkdir -p ${out_d} ; fi

################################################
# Each pop
################################################

for pop in ${pops[@]} ; do
    echo ${pop} ${Csq}
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_${pop}.phe"

    if [ "${Csq}" == "redo" ] ; then
    # based on the notebook 6c_LD_indep_check, we try to apply LD pruning again
        out="${out_d}/$(basename ${pfile}).${pop}.${Csq}.bool"
        plink_indep_pairwise --keep ${keep} --out ${out} \
        --extract ${out_d}/ukb24983_cal_hla_cnv.${pop}.redo.lst
        mv ${out}.log ${out}.prune.log

    elif [ "${Csq}" == "QCed" ] ; then
    # focus on the QCed variants and apply LD pruning without considering the consequence field
        out="${out_d}/$(basename ${pfile}).${pop}.${Csq}.bool"
        cat_or_zcat ${QCed_vars_f} | awk -v FS='\t' '(NR>1){print $3}' \
        | plink_indep_pairwise --keep ${keep} --out ${out} --extract /dev/stdin
        mv ${out}.log ${out}.prune.log

    elif [ "${Csq}" == "" ] ; then
    # apply LD pruning for the entire dataset
        out="${out_d}/$(basename ${pfile}).${pop}.bool"
        plink_indep_pairwise --keep ${keep} --out ${out}
        mv ${out}.log ${out}.prune.log
        
    else
    # apply LD pruning for the specified consequence type, given the LD pruned 
    # results for more severe consequence types
        out="${out_d}/$(basename ${pfile}).${pop}.${Csq}.bool"
        show_var_lst ${Csq} \
        | plink_indep_pairwise --keep ${keep} --out ${out} --extract /dev/stdin
        mv ${out}.log ${out}.prune.log
    fi
done
