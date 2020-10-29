#!/bin/bash
set -beEuo pipefail

n_batch=900

data_d='/scratch/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020'

{
    zcat ${data_d}/liftUnder/UKBexomeOQFE.001.tsv.gz | egrep '#'
    seq -w ${n_batch} | while read idx_pad ; do
        zcat ${data_d}/liftUnder/UKBexomeOQFE.${idx_pad}.tsv.gz | egrep -v '#'
    done 
} | bgzip -f -l9 -@6 > ${data_d}/UKBexomeOQFE.hg19.tsv.gz

seq -w ${n_batch} | while read idx_pad ; do
    zcat ${data_d}/liftUnder/UKBexomeOQFE.${idx_pad}.unmapped.txt
done | bgzip -f -l9 -@6 > ${data_d}/UKBexomeOQFE.hg19.unmapped.txt.gz

