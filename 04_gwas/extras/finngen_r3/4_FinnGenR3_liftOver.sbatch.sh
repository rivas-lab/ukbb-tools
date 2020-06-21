#!/bin/bash
set -beEuo pipefail

index_file="4_FinnGenR3_liftOver.input.lst"

batch_idx=${SLURM_ARRAY_TASK_ID:=361}
batch_size=5

# ml load snpnet_yt

start_idx=$( perl -e "print(int(((${batch_idx} - 1) * ${batch_size}) + 1))" )
end_idx=$(   perl -e "print(int(  ${batch_idx}      * ${batch_size}))" )

cat ${index_file} | egrep -v '^#' | awk -v s=${start_idx} -v e=${end_idx} '(s <= NR && NR <= e){print $1}' |
while read f ; do
    if [ -f $f ] ; then
        bash 4_FinnGenR3_liftOver.sh $f
    fi
done

exit 0
#######################

# usage
ml load R/3.6 gcc
sbatch -p mrivas --qos=high_p --nodes=1 --mem=8000 --cores=2 --time=5:00:00 --job-name=FGr3 --output=logs/FGr3.%A_%a.out --error=logs/FGr3.%A_%a.err --array=1-361 4_FinnGenR3_liftOver.sbatch.sh
Submitted batch job @@@@@

# file location
Input:   /scratch/groups/mrivas/public_data/summary_stats/finngen_r3/summary_stats
Results: /scratch/groups/mrivas/public_data/summary_stats/finngen_r3/summary_stats_hg19
