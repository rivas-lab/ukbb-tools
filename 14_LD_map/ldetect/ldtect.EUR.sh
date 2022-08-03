#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f "${0}")
SRCDIR=$(dirname "${SRCNAME}")

# assign LDetect LD block idx
# https://bitbucket.org/nygcresearch/ldetect-data

ukb24983_d="/oak/stanford/groups/mrivas/ukbb24983"

ldtect_f="/oak/stanford/groups/mrivas/users/ytanigaw/repos/nygcresearch/ldetect-data/EUR/fourier_ls-all.bed"

pvar_f="${ukb24983_d}/array-combined/pgen/ukb24983_cal_hla_cnv.pvar.gz"

pvar_ldtect_f="${ukb24983_d}/array-combined/ldetect/ukb24983_cal_hla_cnv.ldtect.EUR.tsv.gz"

{

    tabix -H ${pvar_f} \
    | awk -v OFS='\t' -v ld_block_idx="ld_block_idx" -v ld_block_name="ld_block_name" \
        '{print $0, ld_block_idx, ld_block_name}'

    less "${ldtect_f}" | awk '(NR>1){print NR-1, $1, $2, $3}' | sed -e 's/chr//g' \
    | while read -r ld_block_idx chrom pos_s pos_e ; do
        tabix ${pvar_f} ${chrom}:${pos_s}-${pos_e} \
        | awk -v OFS='\t' -v ld_block_idx=${ld_block_idx} -v ld_block_name=${chrom}:${pos_s}-${pos_e} \
        '{print $0, ld_block_idx, ld_block_name}'
    done
} | bgzip -l9 -@6 > ${pvar_ldtect_f}
