#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

ml load R/3.6 gcc

data_d="/oak/stanford/groups/mrivas/ukbb24983/exome/annotation/20201025_exome_oqfe_2020"

prev_f="${data_d}/ukb24983_exomeOQFE.annotation.20201217.compact.tsv.gz"
new_f="${data_d}/ukb24983_exomeOQFE.annotation.20210108.compact.tsv.gz"

Rscript ${SRCNAME%.sh}.R ${prev_f} ${new_f%.gz}
bgzip -l9 -@6 ${new_f%.gz}

prev_f="${data_d}/ukb24983_exomeOQFE.annotation.20201217.tsv.gz"
new_f="${data_d}/ukb24983_exomeOQFE.annotation.20210108.tsv.gz"

Rscript ${SRCNAME%.sh}.R ${prev_f} ${new_f%.gz}
bgzip -l9 -@6 ${new_f%.gz}
