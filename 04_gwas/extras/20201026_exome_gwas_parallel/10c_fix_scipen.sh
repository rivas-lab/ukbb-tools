#!/bin/bash
set -beEuo pipefail

pop='metal'
chrom=$1

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"
gz_f=${data_d}/ukb24983_exomeOQFE.${pop}.chr${chrom}.p1e-3.tsv.gz

ml load R/3.6 gcc

if [ ! -f ${gz_f%.gz} ] ; then
    Rscript 10c_fix_scipen.R ${gz_f} ${gz_f%.gz}
fi

bgzip -@6 -f ${gz_f%.gz}

if [ ! -f ${gz_f}.tbi ] ; then
    tabix -s1 -b2 -e2 -c'#' ${gz_f}
fi

exit 0
########################
for chr in X 1 2 5 8 9 10 11 12 17 19 ; do
for pop in 'metal' ; do
    sbatch -p mrivas --qos=high_p --nodes=1 --mem=32000 --cores=2 --time=2:00:00 \
    --job-name=fix_scipen.${chr} --output=logs/fix_scipen.${chr}.%A.out --error=logs/fix_scipen.${chr}.%A.err \
    10c_fix_scipen.sh ${chr}
done
done

