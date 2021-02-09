#!/bin/bash
set -beEuo pipefail

ml load rclone

src_d="/oak/stanford/groups/mrivas/ukbb24983/exome/gwas_50k/freeze/20201106/"
dst_d=$(echo ${src_d} | sed -e "s%/oak/stanford/groups/mrivas/%gdrive:rivas-lab/%g")

find ${src_d} -type f | sort | egrep 'related|s_asian|others|white_british' | while read f ; do 
    echo $f
    rclone copy $f ${dst_d}/
done

exit 0
