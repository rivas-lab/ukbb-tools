#!/bin/bash
set -beEuo pipefail

pop=$1

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"
out_f=${data_d}/ukb24983_exomeOQFE.${pop}.p1e-3.tsv.gz

if [ ! -f ${out_f} ] ; then
    {
        ! zcat ${data_d}/ukb24983_exomeOQFE.${pop}.chr1.p1e-3.tsv.gz | head -n1

        for chrom in $(seq 1 22) X Y ; do
            zcat ${data_d}/ukb24983_exomeOQFE.${pop}.chr${chrom}.p1e-3.tsv.gz \
            | egrep -v '^#'
        done
    } | bgzip -l9 -@6 > ${out_f}
fi

if [ ! -f ${out_f}.tbi ] ; then
    tabix -s1 -b2 -e2 -c'#' ${out_f}
fi
