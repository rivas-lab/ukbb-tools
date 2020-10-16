#!/bin/bash
set -beEuo pipefail

ml load plink2


pfile="/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv"
QCed_vars_f="/oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/ld_indep_20201015/ukb24983_cal_hla_cnv.var_QCed.tsv.gz"
out_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/annotation/ld_indep_20201015/plink_output"
pops=('white_british')

plink_wrapper () {
    plink2 --silent --threads 6 --memory 60000 --pfile ${pfile} vzs $@
}

plink_indep_pairwise () {
    plink_wrapper --allow-extra-chr --indep-pairwise 50 5 0.5 $@
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

if [ ! -d ${out_d} ] ; then mkdir -p ${out_d} ; fi

################################################
# Each pop
################################################

for pop in ${pops[@]} ; do
    echo ${pop}
    keep="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_${pop}.phe"
    Csq="ptv"


    out="${out_d}/$(basename ${pfile}).${pop}.${Csq}.bool"

    cat_or_zcat ${QCed_vars_f} | awk -v FS='\t' -v Csq=${Csq} '($NF == Csq){print $3}' \
    | plink_indep_pairwise --keep ${keep} --out ${out} --extract /dev/stdin
    mv ${out}.log ${out}.prune.log
done
