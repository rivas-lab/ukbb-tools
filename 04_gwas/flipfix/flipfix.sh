#!/bin/bash
set -beEuo pipefail

# Given a GWAS summary statistics in PLINK format, 
# check the "REF" col against the reference sequence in a FASTA file 

. "$(dirname $(readlink -f $0))/flip_misc.sh"

in_sumstats=$(readlink -f $1)
if [ $# -gt 1 ] ; then ref_fa=$(readlink -f $2) ; else ref_fa="${REF_HG19_FA}" ; fi

tmp_dir_root=/dev/shm/u/$USER
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

#check_flip $in_sumstats $ref_fa ${tmp_dir}
fix_flip ${in_sumstats} ${ref_fa} ${tmp_dir}

