#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
   
pvar="/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE.pvar.gz"
bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/09_liftOver/flipcheck.sh \
    ${pvar} \
    | bgzip -l9 -@6 > ${pvar%.pvar.gz}.flipcheck.hg19.tsv.gz

