#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc

data_d="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc_upset"

! tabix -l ${data_d}/ukb37855_ukb40831_icd.annot.tsv.gz | grep -v NA | parallel --eta -j6 'bash 3_upset_plot.sh {}'

for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' ; do

echo "=============================================================="
echo "${pop}"
echo "=============================================================="
! tabix -l ${data_d}/ukb37855_ukb40831_icd.annot.tsv.gz | grep -v NA | parallel --eta -j6 "bash 3_upset_plot.sh {} ${pop}"

done

find ${data_d} -size 3611c -name "*.pdf" -exec rm {} \;
