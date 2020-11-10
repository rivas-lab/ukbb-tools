#!/bin/bash
set -beEuo pipefail

data_d=/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE

cd ${data_d}

find 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' \
    -name "*.glm.logistic.hybrid.gz" -o -name "*.glm.linear.gz" \
    | sort -V \
    | awk -v d=${data_d} '{ print d "/" $0 }'
