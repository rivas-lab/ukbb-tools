#!/bin/bash
set -beEuo pipefail

ml load plink/1.90b6.21

out="/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE"

zcat /scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen/ukb24983_cal_hla_cnv_exomeOQFE.pvar.gz \
| awk '(NR>1){print $3}' \
| plink --memory 60000 --threads 8 \
--bfile /scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE \
--extract /dev/stdin \
--bmerge /scratch/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv \
--keep-allele-order \
--mac 1 \
--make-bed \
--out ${out}

mv ${out}.log ${out}.bed.log


exit 0

## job submission instruction
sbatch -p mrivas --qos=high_p --time=1-0:0:00 --mem=64000 --nodes=1 --cores=8 --job-name=merge --output=logs/merge.%A.out --error=logs/merge.%A.err 3_merge.sh
