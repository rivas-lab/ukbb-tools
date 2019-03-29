#!/bin/bash 
set -beEu

out_dir=$(readlink -f $1)

sed_key="$(echo  ${out_dir} | sed -e "s%rg%munged%g")/"

echo "#p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se" | tr " " "\t"
find ${out_dir} -type f | sort | while read f ; do 
    cat $f | tail -n4 | head -n1 | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'
done | sed -e "s%${sed_key}%%g" \
| sed -e 's/ukb24983_v2.//g' \
| sed -e 's/.sumstats.gz//g' \
| sed -e 's/.genotyped//g'

