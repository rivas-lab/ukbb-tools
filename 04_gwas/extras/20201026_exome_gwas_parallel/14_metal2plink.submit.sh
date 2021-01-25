#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

src="/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/18_metal/metal_to_plink.sh"

jn="metal2plink"
job_idx="14_metal2plink.submit.job.lst"

find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/metal/ -name "*.tsv.gz" | sort -V > ${job_idx}

sbatch -p mrivas,normal --time=1:0:00 --mem=6000 --nodes=1 --cores=1 \
--job-name=${jn} --output=logs/${jn}.%A_%a.out --error=logs/${jn}.%A_%a.err \
--array=1-339 ${parallel_sbatch_sh} ${src} ${job_idx} 10
