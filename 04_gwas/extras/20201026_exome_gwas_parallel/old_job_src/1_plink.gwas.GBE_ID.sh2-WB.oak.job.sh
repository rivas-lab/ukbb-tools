#!/bin/bash
set -beEuo pipefail

GBE_ID=$1
if [ $# -gt 1 ] ; then batch_idx_s=$2 ; else batch_idx_s=11 ; fi

for pop in 'white_british' ; do

out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch

for batch_idx in $(seq ${batch_idx_s} 100) ; do

echo $batch_idx

! bash 1_plink.gwas.sh \
--GBE_ID ${GBE_ID} \
--pop ${pop} \
--pfile /scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
${batch_idx} \
${out_d}

done
done

