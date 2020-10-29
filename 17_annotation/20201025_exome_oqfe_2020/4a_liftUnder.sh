#!/bin/bash
set -beEuo pipefail

idx=$1
idx_pad=$(perl -e "print(sprintf('%03d', ${idx}))")

# data_d='/scratch/groups/mrivas/ukbb24983/exome-oqfe2020-annotation'
data_d='/scratch/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020'

bash /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/09_liftOver/liftOver_wrapper.sh \
    --threads 1 \
    --src_genome hg38 --dst_genome hg19 \
    ${data_d}/input/split.UKBexomeOQFE.pvar.body.${idx_pad}.pvar \
    ${data_d}/liftUnder/UKBexomeOQFE.${idx_pad}.tsv \
    ${data_d}/liftUnder/UKBexomeOQFE.${idx_pad}.unmapped.txt

