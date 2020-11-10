#!/bin/bash
set -beEuo pipefail

pop=$1
GBE_ID=$2
if [ $# -gt 2 ] ; then batch_idx_s=$3 ; else batch_idx_s=1 ; fi
if [ $# -gt 3 ] ; then batch_idx_e=$4 ; else batch_idx_e=100 ; fi
if [ $# -gt 4 ] && [ $5 == "oak" ] ; then oak="TRUE" ; else oak="FALSE" ; fi

if [ "${oak}" == "TRUE" ] ; then
    out_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE_scg/${pop}-batch
else
    out_d=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/${pop}-batch
fi

for batch_idx in $(seq ${batch_idx_s} ${batch_idx_e}) ; do

echo $batch_idx

! bash 1_plink.gwas.sh \
--GBE_ID ${GBE_ID} \
--pop ${pop} \
--pfile /scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
${batch_idx} \
${out_d}

done

exit 0
###############################################

bash 1_plink.gwas.wrapper.sh white_british FH1044

cat 3a_check_results.BINs.20201102-091358.tsv \
| egrep -v '#' | awk '($3==99){print $1, $2}' \
| while read pop GBE_ID ; do 
    echo ${GBE_ID} ${pop}
    sbatch -p mrivas --qos=high_p --nodes=1 \
    --mem=8000 --cores=2 --time=6:00:00 \
    --job-name=gwas.${GBE_ID} \
    --output=logs_scratch/gwas.${GBE_ID}.%A.out \
    --error=logs_scratch/gwas.${GBE_ID}.%A.err \
    1_plink.gwas.wrapper.sh ${pop} ${GBE_ID}
done
