#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)

ml load python/2.7 R/3.6 gcc

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
############################################################

idx=$1
idx_pad=$(perl -e "print(sprintf('%04d', ${idx}))")

assembly=GRCh37

split_pvar=/scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003/input_20201006/split.cnv.20201006.unfinished.pvar.body.${idx_pad}.pvar
split_seq_pvar=${split_pvar%.pvar}.seq.pvar
vep_in_vcf=${split_pvar%.pvar}.vcf
vep_out=$(dirname $(dirname ${vep_in_vcf}))/output_vep_20201006/$(basename ${vep_in_vcf%} .vcf).vep101-loftee

# reference data
public_d="/oak/stanford/groups/mrivas/public_data"
vep_data="${public_d}/vep/20200912"
loftee_data=$(dirname ${vep_data})/20201002_loftee_data


if [ ! -d $(dirname ${vep_out}) ] ; then mkdir -p $(dirname ${vep_out}) ; fi


if [ ! -s ${vep_out}.tsv ] ; then
    # Allele strings are too long in CNV dataset. Let's clean them up.    
    Rscript simplify_cnv_tsv.R ${vep_out}.tmp.tsv ${vep_out}.tsv 

    rm ${vep_out}.tmp.tsv
    bgzip ${vep_out}.vcf
fi
