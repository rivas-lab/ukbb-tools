#!/bin/bash
set -beEuo pipefail

source /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/14_LD_map/array-exome-combined/20210112/0_parameters.sh

ml load plink2 R/3.6 gcc

# Csq="ptv"
if [ $# -gt 0 ] ; then Csq=$1 ; else Csq="" ; fi

pfile="/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE"
pfile="/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen_20201217_bug/ukb24983_cal_hla_cnv_exomeOQFE"

QCed_vars_f="${ldmap_d}/ukb24983_cal_hla_cnv_exomeOQFE.input.variants.tsv.gz"
out_d="${ldmap_d}/plink_output"
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
    local Csq=$1
    local rsrc="2c_LD_indep_var_filter.R"
    
    if   [ ${Csq} == "ptv" ] ; then
        Rscript ${rsrc} ${Csq} 'none' /dev/stdout
    elif [ ${Csq} == "pav" ] ; then
        Rscript ${rsrc} ${Csq} 'ptv' /dev/stdout
    elif [ ${Csq} == "pcv" ] ; then
        Rscript ${rsrc} ${Csq} 'ptv,pav' /dev/stdout
    elif [ ${Csq} == "utr" ] ; then
        Rscript ${rsrc} ${Csq} 'ptv,pav,pcv' /dev/stdout
    elif [ ${Csq} == "intron" ] ; then
        Rscript ${rsrc} ${Csq} 'ptv,pav,pcv,utr' /dev/stdout
    elif [ ${Csq} == "others" ] ; then
        Rscript ${rsrc} ${Csq} 'ptv,pav,pcv,utr,intron' /dev/stdout
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
    
    keep_f=$(echo ${keep_f_template} | sed -e "s/__POPULATION__/${pop}/g")

    if [ "${Csq}" == "redo" ] ; then
    # based on the notebook 6c_LD_indep_check, we try to apply LD pruning again
        out="${out_d}/$(basename ${pfile}).${pop}.${Csq}.bool"
        plink_indep_pairwise --keep ${keep_f} --out ${out} \
        --extract ${out_d}/ukb24983_cal_hla_cnv.${pop}.redo.lst
        mv ${out}.log ${out}.prune.log

    elif [ "${Csq}" == "QCed" ] ; then
    # focus on the QCed variants and apply LD pruning without considering the consequence field
        out="${out_d}/$(basename ${pfile}).${pop}.${Csq}.bool"
        cat_or_zcat ${QCed_vars_f} | awk -v FS='\t' '(NR>1){print $3}' \
        | plink_indep_pairwise --keep ${keep_f} --out ${out} --extract /dev/stdin
        mv ${out}.log ${out}.prune.log

    elif [ "${Csq}" == "" ] ; then
    # apply LD pruning for the entire dataset
        out="${out_d}/$(basename ${pfile}).${pop}.bool"
        plink_indep_pairwise --keep ${keep_f} --out ${out}
        mv ${out}.log ${out}.prune.log
        
    else
    # apply LD pruning for the specified consequence type, given the LD pruned 
    # results for more severe consequence types
        out="${out_d}/$(basename ${pfile}).${pop}.${Csq}.bool"
        show_var_lst ${Csq} \
        | plink_indep_pairwise --keep ${keep_f} --out ${out} --extract /dev/stdin
        mv ${out}.log ${out}.prune.log
    fi
done
