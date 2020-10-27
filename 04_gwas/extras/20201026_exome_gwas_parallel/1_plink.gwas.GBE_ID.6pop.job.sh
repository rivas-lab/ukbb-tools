#!/bin/bash
set -beEuo pipefail

GBE_ID=$1

for pop in 'e_asian' 'african' 's_asian' 'non_british_white' 'others' 'related' ; do

out_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch

for batch_idx in $(seq 1 100) ; do

echo $batch_idx

! bash 1_plink.gwas.sh \
--GBE_ID ${GBE_ID} \
--pop ${pop} \
--n_batch 100 \
--pfile /scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
${batch_idx} \
${out_d}

done
done

