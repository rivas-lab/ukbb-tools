#!/bin/bash
set -beEuo pipefail

cd /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE
find 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' 'metal' -type f -name "*linear.gz" -o -name "*logistic.hybrid.gz" -o -name "*metal.tsv.gz" | awk -v FS='/' '{print $1}' | sort | uniq -c | sort -k1,1nr | awk -v OFS='\t' '{print $1, $2}'

