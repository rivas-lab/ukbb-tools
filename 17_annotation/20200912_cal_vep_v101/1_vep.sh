#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"

public_d="/oak/stanford/groups/mrivas/public_data"
vep_data="${public_d}/vep/20200912"
data_d=/oak/stanford/groups/mrivas/ukbb24983/cal
pvar=${data_d}/pgen/ukb24983_cal_cALL_v2_hg19.pvar
vep_out=${data_d}/annotation_20200912/ukb24983_cal_cALL_v2_hg19.vep101.noLoF
vep_in_vcf=${vep_out}.input.vcf

if [ ! -d $(dirname ${vep_out}) ] ; then mkdir -p $(dirname ${vep_out}) ; fi

# load the vep module 
# see: https://github.com/rivas-lab/sherlock-modules/tree/master/vep
ml load vep

# convert pvar to vcf
! bash $(dirname ${SRCDIR})/anno_yt/pvar.to.vcf.sh ${pvar} \
| sed -e 's/chrMT/chrM/g' | sed -e 's/chrXY/chrX/g' > ${vep_in_vcf}

# run vep
vep --fork 6 \
--offline --cache \
--allele_number --everything \
--assembly GRCh37 -i ${vep_in_vcf} -o ${vep_out}

# clean-up
mv ${vep_out} ${vep_out}.tsv
rm ${vep_in_vcf}
