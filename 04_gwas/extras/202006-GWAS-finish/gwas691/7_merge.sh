#!/bin/bash
set -beEuo pipefail

map_f=/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas691/6_map.tsv

task_idx=$1

f1=$(cat $map_f | egrep -v '^#' | awk -v nr=${task_idx} '(NR==nr){print $4}')
f2=$(cat $map_f | egrep -v '^#' | awk -v nr=${task_idx} '(NR==nr){print $5}')

bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas369/4_combine.sh $f1 $f2

exit 0
####
# Instructions on the batch job
seq $(cat 6_map.tsv | egrep -v '^#' | wc -l ) > 7_merge.job.lst

sbatch -p mrivas,normal --nodes=1 --mem=4000 --cores=1 --time=1:00:00 --job-name=g691 --output=logs/g691.%A_%a.out --error=logs/g691.%A_%a.err --array=1-180 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 7_merge.sh 7_merge.job.lst 1
