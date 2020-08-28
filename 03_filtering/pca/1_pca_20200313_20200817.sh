#!/bin/bash
set -beEuo pipefail

sqc_d="/oak/stanford/groups/mrivas/ukbb24983/sqc"
sqc_v="20200313"
pca_v="20200817_v3"

# v1 --indep-pairwise 50 5 .5
# v2 --indep-pairwise 50 5 .2
# v3 we use UKB's SNP QC file and focus on the 147604 variats used in their PCA analysis

pop=$1
shift 1

bash 1_pca.sh \
--keep ${sqc_d}/population_stratification_w24983_${sqc_v}/ukb24983_${pop}.phe \
${sqc_d}/population_stratification_w24983_${sqc_v}/pca_${pca_v}/ukb24983_${pop} \
$@
