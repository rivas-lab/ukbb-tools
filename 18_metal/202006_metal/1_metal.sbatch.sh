#!/bin/bash
set -beEuo pipefail

info_file="../../05_gbe/phenotype_info.tsv"

batch_idx=${SLURM_ARRAY_TASK_ID:=1}
if [ $# -gt 0 ] ; then batch_idx=$1 ; fi
batch_size=10

# ml load snpnet_yt

start_idx=$( perl -e "print(int(((${batch_idx} - 1) * ${batch_size}) + 1))" )
end_idx=$(   perl -e "print(int(  ${batch_idx}      * ${batch_size}))" )

cat ${info_file} | egrep -v '^#' | awk -v s=${start_idx} -v e=${end_idx} '(s <= NR && NR <= e){print $1}' |
while read GBE_ID ; do
    echo $GBE_ID >&2
    bash 1_metal.sh $GBE_ID
done

exit 0
#######################

# usage
ml load snpnet_yt
sbatch -p mrivas --qos=high_p --nodes=1 --mem=4000 --cores=1 --time=3:00:00 --job-name=metal --output=logs/metal.%A_%a.out --error=logs/metal.%A_%a.err --array=1-402 1_metal.sbatch.sh
Submitted batch job 2499937

# file location
Results: /oak/stanford/groups/mrivas/ukbb24983/array-combined/metal/20200616/
