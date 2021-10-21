#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

ukb_d='/oak/stanford/groups/mrivas/ukbb24983'

sqc_files=(
pca1_refinement
pca2_covars
PCA.plots.tar.gz
ukb24983_white_british.phe
ukb24983_non_british_white.phe
ukb24983_african.phe
ukb24983_s_asian.phe
ukb24983_e_asian.phe
ukb24983_others.phe
ukb24983_related.phe
)

cd ${ukb_d}/sqc/population_stratification_w24983_20211020
for f in ${sqc_files[@]} ; do
    set -x
    echo ln -s ../population_stratification_w24983_20200828/${f} .
    set +x
done

cd ${ukb_d}/phenotypedata/master_phe

set -x

ln -s master.20210129.phe.info.tsv master.20211020.phe.info.tsv

ln -sf master.20211020.phe.info.tsv master.phe.info.tsv
ln -sf master.20211020.phe          master.phe

zcat master.20211020.phe.gz > master.20211020.phe
zstd -9 master.20211020.phe
