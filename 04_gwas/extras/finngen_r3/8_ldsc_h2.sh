#!/bin/bash
set -beEuo pipefail

FinnGen=$1

finngen_d="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3"
FinnGen_f="${finngen_d}/summary_stats_hg19_ldsc_munge/finngen_r3_${FinnGen}.hg19.sumstats.gz"
out_f="${finngen_d}/ldsc/h2/${FinnGen}"

src="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_h2.sh"

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

if [ ! -f ${out_f}.log ] ; then
    bash ${src} --scratch ${FinnGen_f} ${out_f}
fi
