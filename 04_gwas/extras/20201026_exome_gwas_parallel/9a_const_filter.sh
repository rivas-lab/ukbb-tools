#!/bin/bash
set -beEuo pipefail

in_f=$1

# in_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/white_british/ukb24983_exomeOQFE.INI25521.glm.linear.gz

nCores=1
if [ $# -gt 1 ] ; then nCores=$2 ; fi

out_f=$(dirname $(dirname $(dirname ${in_f})))/$(basename $(dirname $(dirname ${in_f})))_filtered/$(basename $(dirname ${in_f}))/$(basename ${in_f})

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

###

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

if [ ! -s ${out_f} ] ; then

cat_or_zcat ${in_f} \
| awk -v FS='\t' -v const="CONST_OMITTED_ALLELE" '$NF != const' \
| bgzip -l9 -@ ${nCores} > ${out_f}

fi

# echo $out_f
exit 0
#######################
find /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE -type f \
    -name "*glm.logistic.hybrid.gz" -o -name "*.glm.linear.gz" \
    | sort -V > 9_filter.$(date +%Y%m%d-%H%M%S).lst

20013
9_filter.sh

sbatch -p mrivas,normal --time=2:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=filter --output=logs/filter.%A_%a.out --error=logs/filter.%A_%a.err \
    --array=1-401 ${parallel_sbatch_no_err_check_sh} 9_filter.sh 9_filter.20201107-151533.lst 50

