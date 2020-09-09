#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc inkscape

data_d="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc_upset"

! tabix -l ${data_d}/ukb37855_ukb40831_icd.annot.tsv.gz | grep -v NA | parallel --eta -j6 'bash 3_upset_plot.sh {}'

! tabix -l ${data_d}/ukb37855_ukb40831_icd.annot.tsv.gz | grep -v NA | parallel --eta -j6 "bash 3_upset_plot.sh {1} {2}" :::: /dev/stdin ::: white_british non_british_white african s_asian e_asian related others

find ${data_d} -size 3611c -name "*.pdf" -exec rm {} \;

