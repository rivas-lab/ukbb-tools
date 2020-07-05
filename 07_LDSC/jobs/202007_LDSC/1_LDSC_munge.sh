#!/bin/bash
set -beEuo pipefail

# ml load ukbb-tools/20200630 R/3.6 gcc

in_f=$1
#in_f=/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/african/ukb24983_v2_hg19.BIN100020.array-combined.glm.logistic.hybrid.gz

out_f=$(dirname $(dirname $(dirname $(dirname $in_f))))/ldsc/$(basename $(dirname $in_f))/$(basename $in_f | sed -e "s/.glm.linear.gz//g" | sed -e "s/.glm.logistic.hybrid.gz//g" )

if [ ! -d $(dirname $out_f) ] ; then mkdir -p $(dirname $out_f) ; fi

if [ ! -f ${out_f}.log ] || [ ! -f ${out_f}.sumstats.gz ] ; then
    bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch $in_f $out_f
fi

exit 0
###########################
# usage

ml load R/3.6 gcc

find /oak/stanford/groups/mrivas/users/ytanigaw/data/.20200605_NGH_UKB_phase1_metabolomics/gwas/ -name "*.gz" | sort | tee 1_LDSC_munge_input.lst

sbatch -p mrivas,owners,normal --nodes=1 --mem=8000 --cores=1 --time=1:00:00 --job-name=ldsc_munge --output=logs/ldsc_munge.%A_%a.out --error=logs/ldsc_munge.%A_%a.err --array=1-996 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 1_LDSC_munge.sh 1_LDSC_munge_input.lst 2

