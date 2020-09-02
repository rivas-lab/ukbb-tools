#!/bin/bash
set -beEuo pipefail

sqc_d="/oak/stanford/groups/mrivas/ukbb24983/sqc"
sqc_v="20200828"
pca_v="pca2_covars"

for pop in 's_asian' ; do

bash ../pca/1_pca.sh \
--keep ${sqc_d}/population_stratification_w24983_${sqc_v}/${pca_v}/ukb24983_${pop}.phe \
${sqc_d}/population_stratification_w24983_${sqc_v}/${pca_v}/ukb24983_${pop}

done
