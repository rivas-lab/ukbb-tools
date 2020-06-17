#!/bin/bash
set -beEuo pipefail

in_file=$1
out_file=$(dirname $(dirname ${in_file}))/summary_stats_hg19/$(basename ${in_file} .gz).hg19

if [ ! -f ${out_file}.gz ] ; then

    Rscript 4_FinnGenR3_liftOver.R ${in_file} ${out_file}
    bgzip -l9 -@2 ${out_file}
    tabix -c '#' -s 1 -b2 -e2 ${out_file}.gz

    echo ${out_file}.gz

fi
