#!/bin/bash
set -beEuo pipefail

data_d='/scratch/groups/mrivas/ukbb24983/cal/ldmap/ldmap_20201018'
pop='white_british'
r2='0.5'
ld_map_f="${data_d}/ukb24983_cal.${pop}.ld_map.tsv.gz"

zcat ${ld_map_f} | awk -v r2=${r2} 'NR==1 || $NF >= r2' | bgzip -l9 -@6 > ${ld_map_f%.tsv.gz}.${r2}r2.tsv.gz

tabix -c '#' -s 1 -b 2 -e 5 ${ld_map_f%.tsv.gz}.${r2}r2.tsv.gz

