#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

pvcf=/scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/download/ukb23156_cY_b0_v1.vcf.gz

zcat ${pvcf} | egrep '^##' | egrep -v '^##contig' > pvcf.header.txt
