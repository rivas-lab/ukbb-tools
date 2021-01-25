#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

source /oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/14_LD_map/array-exome-combined/20210112/0_parameters.sh

in_f="${ldmap_d}/plink_output/ukb24983_cal_hla_cnv_exomeOQFE.__POP__.__CSQ__.bool.prune.__EXT__"
out_f="${ldmap_d}/ukb24983_cal_hla_exomeOQFE.__POP__.bool.__EXT__"
csqs=(ptv pav pcv utr intron others)

pop="white_british"
out_in=$( echo ${out_f} | sed -e "s/__POP__/${pop}/g" | sed -e "s/__EXT__/in/g" )
out_log=$(echo ${out_f} | sed -e "s/__POP__/${pop}/g" | sed -e "s/__EXT__/log/g" )

# LD independent set of variants in the array-combined dataset
zcat ${ldmap_d}/ukb24983_cal_hla_cnv_exomeOQFE.input.variants.tsv.gz \
    | awk -v FS='\t' '($8 == "TRUE"){print $3}' >> ${out_in}

for csq in ${csqs[@]} ; do
    in=$( echo ${in_f} | sed -e "s/__POP__/${pop}/g" | sed -e "s/__CSQ__/${csq}/g" | sed -e "s/__EXT__/in/g" )
    log=$(echo ${in_f} | sed -e "s/__POP__/${pop}/g" | sed -e "s/__CSQ__/${csq}/g" | sed -e "s/__EXT__/log/g" )
    cat ${in}        >> ${out_in}
    echo "## ${log}" >> ${out_log}
    cat ${log}       >> ${out_log}
done

bgzip -l9 -@6 ${out_in}
bgzip -l9 -@6 ${out_log}
