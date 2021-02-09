#!/bin/bash
set -beEuo pipefail

ml load rclone

src_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/freeze/20201106"
dst_d=$(echo ${src_d} | sed -e "s%/oak/stanford/groups/mrivas/%gdrive:rivas-lab/%g")

find ${src_d} -type f | grep metal | grep -v tar | sort | while read f ; do 
    echo $f
    rclone copy $f ${dst_d}/
done

exit 0
