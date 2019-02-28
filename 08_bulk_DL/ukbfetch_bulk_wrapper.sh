#!/bin/bash
set -beEuo pipefail

bulk_file=$1
out_dir=$(readlink -f $2)

n_lines=$(cat ${bulk_file} | wc -l)
n_batches=$( perl -e "print(int((${n_lines} + 999)/1000))" )

parallel -j 4 $(dirname $(readlink -f $0))/$(basename $0 .sh)_sub.sh ${bulk_file} {} ${out_dir} :::  $(seq 1 ${n_batches}) 

