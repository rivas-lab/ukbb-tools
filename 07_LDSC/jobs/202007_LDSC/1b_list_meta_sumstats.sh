#!/bin/bash
set -beEuo pipefail

metal_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/metal/20200717"
out_f="1_LDSC_munge.$(date +%Y%m%d-%H%M%S).metal.job.lst"

find ${metal_d} -type f -name "*.metal.tsv.gz"| sort > ${out_f}
echo ${out_f}
