#!/bin/bash
set -beEuo pipefail

ldmap_d="/scratch/groups/mrivas/ukbb24983/array-exome-combined/ldmap/ldmap_20201223"
pop="white_british"
out_prefix="${ldmap_d}/ukb24983_cal_hla_exomeOQFE.${pop}"
out_tsvgz="${out_prefix}.ld_map.tsv.gz"

zcat ${ldmap_d}/ukb24983_cal_hla_exomeOQFE.${pop}.chrMT.ld_map.tsv.gz | egrep '^#' > ${out_tsvgz%.gz}

for chr_idx in $(seq 1 26) ; do
    chr=$(echo ${chr_idx} | sed -e 's/23/X/g' | sed -e 's/24/Y/g' | sed -e 's/25/XY/g' | sed -e 's/26/MT/g')
    chr_prefix="${ldmap_d}/ukb24983_cal_hla_exomeOQFE.${pop}.chr${chr}"
    
    zcat ${chr_prefix}.ld_map.tsv.gz | egrep -v '^#' >> ${out_tsvgz%.gz}
    {
        echo "## ${chr_prefix}.ld_map.log ##"
        cat "${chr_prefix}.ld_map.log"
        echo ""    
    } >> ${out_tsvgz%.tsv.gz}.log
done

bgzip -l9 -@6 ${out_tsvgz%.gz}
tabix -c '#' -s 1 -b 2 -e 5 ${out_tsvgz}
