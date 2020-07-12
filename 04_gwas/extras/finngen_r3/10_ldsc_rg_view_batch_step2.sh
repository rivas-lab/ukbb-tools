#!/bin/bash
set -beEuo pipefail

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
n_batchs=1000
if [ $# -gt 0 ] ; then batch_idx=$1 ; fi
if [ $# -gt 1 ] ; then n_batchs=$2 ; fi

log_dir="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/ldsc/UKB_WB_rg"
out_dir="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/ldsc/UKB_WB_rg.20200711-171504"

ml load ukbb-tools

in_files="${out_dir}/LDSC.rg.lst"
batch_out="${out_dir}/LDSC.rg.batch${batch_idx}.tsv"

n_in_files=$(cat ${in_files} | wc -l)
n_in_batch=$(perl -e "print(int(0.99999 + ${n_in_files}/${n_batchs}))")
start_idx=$(perl -e "print(${n_in_batch} * (${batch_idx}-1) + 1)")
end_idx=$(perl -e "print(${n_in_batch} * ${batch_idx})")

echo ${start_idx} ${end_idx} ${in_files} >&2

! cat ${in_files} | egrep -v '#' | awk -v s=${start_idx} -v e=${end_idx} 's <= NR && NR <= e' | ldsc_rg_view.sh -l /dev/stdin \
| sed -e "s%/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/summary_stats_hg19_ldsc_munge/finngen_r3_%%g" \
| sed -e "s/.hg19.sumstats.gz//g" \
| sed -e "s%/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/white_british/ukb24983_v2_hg19.%%g" \
| sed -e "s/.array-combined.sumstats.gz//g" \
| sed -e 's/^p1/#p1/g' > ${batch_out}

echo ${batch_out}

exit 0

########################
# instruction
########################
ml load resbatch

seq 1000 > 1.1000.lst

sbatch -p mrivas --qos=high_p --time=30:00 --mem=6000 --nodes=1 --cores=1 --job-name=rg_agg --output=logs/rg_agg.%A_%a.out --error=logs/rg_agg.%A_%a.err --array=1-1000 $parallel_sbatch_sh 10_ldsc_rg_view_batch_step2.sh 1.1000.lst 1
