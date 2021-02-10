#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

ml load R/3.6 gcc

full_f=$1
p_thr='1e-3'
# helper script
gwas_filter_by_p_R='/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/gwas_filter_by_p.R'

filtered_f=$(dirname $(dirname ${full_f}))/$(basename $(dirname ${full_f}))_p${p_thr}/$(basename ${full_f})

if [ ! -d $(dirname ${filtered_f}) ] ; then mkdir -p $(dirname ${filtered_f}) ; fi

if [ ! -f ${filtered_f} ] ; then
    if [ ! -f ${filtered_f%.gz} ] ; then
        Rscript ${gwas_filter_by_p_R} ${full_f} ${filtered_f%.gz} ${p_thr}
    fi
    bgzip -l9 ${filtered_f%.gz}
fi

exit 0

#############################
# instruction for job-submission on Sherlock

find /scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE \
-type f -name "*.gz" | sort -V > 12c_filter_by_pval.input.lst

sbatch -p mrivas,normal --nodes=1 --mem=8000 --cores=1 --time=12:00:00 --job-name=filter --output=logs/filter.%A_%a.out --error=logs/filter.%A_%a.err --array=1-468 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 12c_filter_by_pval.sh 12c_filter_by_pval.input.lst 50


