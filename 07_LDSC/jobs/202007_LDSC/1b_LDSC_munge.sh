#!/bin/bash
set -beEuo pipefail

ml load R/3.6 gcc

in_f=$1
#in_f=/oak/stanford/groups/mrivas/ukbb24983/array-combined/metal/20200717/INI25466.metal.tsv.gz
GBE_ID=$(basename $in_f .metal.tsv.gz)
pop="metal"
geno_dataset="array-combined"
ldsc_dir="/oak/stanford/groups/mrivas/ukbb24983/${geno_dataset}/ldsc"
out_f="${ldsc_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.${geno_dataset}.sumstats.gz"

if [ ! -d $(dirname $out_f) ] ; then mkdir -p $(dirname $out_f) ; fi

if [ ! -f ${out_f}.log ] || [ ! -f ${out_f}.sumstats.gz ] ; then
    bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch $in_f $out_f
fi

exit 0
###########################
# usage

https://github.com/rivas-lab/ukbb-tools/issues/26

ml load R/3.6 gcc resbatch

sbatch -p mrivas,normal,owners --time=1:00:00 --mem=8000 --nodes=1 --cores=1 --job-name=munge_meta --output=logs/munge_meta.%A_%a.out --error=logs/munge_meta.%A_%a.err --array=1-949 $parallel_sbatch_sh 1b_LDSC_munge.sh 1_LDSC_munge.20200718-134522.metal.job.lst 4
