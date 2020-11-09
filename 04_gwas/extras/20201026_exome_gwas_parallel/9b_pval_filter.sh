#!/bin/bash
set -beEuo pipefail

in_f=$1
p_thr="1e-3"

# in_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_filtered/white_british/ukb24983_exomeOQFE.INI25521.glm.linear.gz

nCores=1
if [ $# -gt 1 ] ; then nCores=$2 ; fi

out_f=$(dirname $(dirname $(dirname ${in_f})))/$(basename $(dirname $(dirname ${in_f})))_p${p_thr}/$(basename $(dirname ${in_f}))/$(basename ${in_f})

###

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

show_header () {
    local file=$1
    ! cat_or_zcat $file | head -n1
    # cat_or_zcat $file | egrep '^#'
}

get_col_idx () {
    local file=$1
    local key=$2
    show_header $file | sed -e "s/^#//g" | tr "\t" "\n" | awk -v key=$key '($0 == key){print NR}'
}

###

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

if [ ! -s ${out_f} ] ; then

p_col=$(get_col_idx ${in_f} "P")

cat_or_zcat ${in_f} \
| awk -v FS='\t' -v p_thr=${p_thr} -v p_col=${p_col} 'NR == 1 || $p_col <= p_thr' \
| bgzip -l9 -@ ${nCores} > ${out_f}

fi

# echo $out_f
exit 0
#######################
find /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_filtered -type f \
    -name "*glm.logistic.hybrid.gz" -o -name "*.glm.linear.gz" \
    | sort -V > 9b_pval_filter.$(date +%Y%m%d-%H%M%S).lst

20013 files

sbatch -p mrivas,normal --time=2:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=filter --output=logs/filter.%A_%a.out --error=logs/filter.%A_%a.err \
    --array=1-401 ${parallel_sbatch_no_err_check_sh} 9b_pval_filter.sh 9b_pval_filter.20201107-222259.lst 50

# metal

find -L /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_filtered/metal -name "*.metal.tsv.gz" \
    | sort -V > 9b_pval_filter.metal.$(date +%Y%m%d-%H%M%S).lst

3381 files

sbatch -p mrivas,normal --time=2:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=filter --output=logs/filter.%A_%a.out --error=logs/filter.%A_%a.err \
    --array=1-846 ${parallel_sbatch_no_err_check_sh} 9b_pval_filter.sh 9b_pval_filter.metal.20201108-093633.lst 4

