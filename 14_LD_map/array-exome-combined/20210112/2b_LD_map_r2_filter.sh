#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

source /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/14_LD_map/array-exome-combined/20210112/0_parameters.sh

ldmap_f="${ldmap_d}/ukb24983_cal_hla_exomeOQFE.white_british.ld_map.tsv.gz"

# filter the LD map file and keep the entries where r2 >= .5

zcat ${ldmap_f} | awk 'NR==1 || $NF>=.5' | bgzip -l9 -@6 > ${ldmap_f%.tsv.gz}.0.5r2.tsv.gz
