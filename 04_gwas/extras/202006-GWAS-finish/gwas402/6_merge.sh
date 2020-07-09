#!/bin/bash
set -beEuo pipefail

map_f=/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas402/20200703-gwas-402.5.map.tsv

task_idx=$1

f1=$(cat $map_f | egrep -v '^#' | awk -v nr=${task_idx} '(NR==nr){print $4}')
f2=$(cat $map_f | egrep -v '^#' | awk -v nr=${task_idx} '(NR==nr){print $5}')

bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/202006-GWAS-finish/gwas369/4_combine.sh $f1 $f2

exit 0
####
# Instructions on the batch job

sbatch -p mrivas,normal --nodes=1 --mem=4000 --cores=1 --time=4:00:00 --job-name=g402 --output=logs/g402.%A_%a.out --error=logs/g402.%A_%a.err --array=1-950 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 6_merge.sh 6_merge.job.lst 15
