#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"

out_fig_d=figs

paths_idx=${SRCDIR}/0_file_names.R

pca1_d=$(cat ${paths_idx} | grep 'pop_refinement_pca' | awk -v FS='=' '{print $NF}' | sed -e 's/,$//g' | sed -e 's/\s//g' | sed -e "s/'//g")
pca2_d=$(cat ${paths_idx} | grep 'pop_specific_pca'   | awk -v FS='=' '{print $NF}' | sed -e 's/,$//g' | sed -e 's/\s//g' | sed -e "s/'//g")

for d in ${pca1_d} ${pca2_d} ; do

    if [ ! -d ${out_fig_d}/$(basename $d) ] ; then mkdir -p ${out_fig_d}/$(basename $d) ; fi

#     find ${d} -maxdepth 1 -name "*.PC1.PC2.png"
    find ${d} -maxdepth 1 -name "*.PC1.PC2.png" -exec cp {} ${out_fig_d}/$(basename $d)/ \;

done
