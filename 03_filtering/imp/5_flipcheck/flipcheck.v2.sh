#!/bin/bash
set -beEuo pipefail

pvar_zst="/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen_v2/ukb24983_cal_hla_cnv_imp.pvar.zst"

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT


tmp_pvar=${tmp_dir}/flipcheck.in.tsv

zstdcat ${pvar_zst} \
| awk '$4 != "N" && $4 != "P" && $5 != "+"' > ${tmp_pvar}

bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/09_liftOver/flipcheck.sh ${tmp_pvar} | awk 'toupper($4) != toupper($NF)' > flipcheck.v2.out
