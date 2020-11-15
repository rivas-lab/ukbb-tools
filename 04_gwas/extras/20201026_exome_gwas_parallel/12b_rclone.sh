#!/bin/bash
set -beEuo pipefail

pop=$1

ml load rclone

echo "[$(date +%Y%m%d-%H%M%S)] ${pop} start"

rclone copy \
/scratch/groups/mrivas/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110/ukb24983_exomeOQFE.${pop}.p1e-3.tsv.gz \
gdrive://rivas-lab/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110

# rclone copy \
# /scratch/groups/mrivas/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110/ukb24983_exomeOQFE.${pop}.p1e-3.tsv.gz.tbi \
# gdrive://rivas-lab/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110

# rclone copy \
# /scratch/groups/mrivas/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110/master_phe_20201002_exomeOQFE_20201110.${pop}.tar \
# gdrive://rivas-lab/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110
echo "[$(date +%Y%m%d-%H%M%S)] ${pop} done"

exit 0
##############################

sbatch -p mrivas --qos=high_p \
--nodes=1 --mem=6000 --cores=1 --time=7-0:00:00 \
--job-name=rclone.${pop} --output=logs/rclone.${pop}.%A.out --error=logs/rclone.${pop}.%A.err \
12b_rclone.sh ${pop}

