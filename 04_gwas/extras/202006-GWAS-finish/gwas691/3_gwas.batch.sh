#!/bin/bash
set -beEuo pipefail
wc_l=$1

cat ../gwas-current-gz-wc.20200629-131527.tsv | awk -v FS='\t' -v wc_l=${wc_l} '($3==wc_l){print $1}' | while read symlink ; do
    pop=$(basename $(dirname $symlink))
    GBE_ID=$(basename $symlink | awk -v FS='.' '{print $2}')
    echo $GBE_ID $pop
    bash 3_gwas.sh $GBE_ID $pop
done
