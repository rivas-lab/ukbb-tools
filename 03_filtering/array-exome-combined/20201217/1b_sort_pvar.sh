#!/bin/bash
set -beEuo pipefail

data_d="/scratch/groups/mrivas/ukbb24983/array-exome-combined/pgen"
in_f="${data_d}/merge_list_pvar/ukb24983_cal_hla_cnv_exomeOQFE.unsorted.pvar"
out_f=${in_f%.unsorted.pvar}.pvar

cat ${in_f} | egrep '^#' | cut -f1-6 > ${out_f}
cat ${in_f} | egrep -v '^#' | sort --parallel 6 -k7,7n -k2,2n -k3,3V | cut -f1-6 >> ${out_f}

zstd -9 ${out_f}
bgzip -@6 -l9 ${out_f}
rm ${in_f}

for ext in gz zst ; do
    cp -a ${out_f}.${ext} $(echo ${out_f}.${ext} | sed -e 's%/scratch/%/oak/stanford/%g')
done

