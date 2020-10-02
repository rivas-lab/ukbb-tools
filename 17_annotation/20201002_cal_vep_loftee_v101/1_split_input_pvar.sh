#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)

pvar="/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19.pvar"
data_d="/scratch/groups/mrivas/ukbb24983/cal/annotation_20201002/input"

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] || [ "${file%.bgz}.bgz" == "${file}" ] ; then
        zcat ${file}
    elif [ "${file%.zst}.zst" == "${file}" ] ; then
        zstdcat ${file}
    else
        cat ${file}
    fi
}

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

if [ ! -d ${data_d} ] ; then mkdir -p ${data_d} ; fi

tmp_pvar_head=${tmp_dir}/original.$(basename ${pvar%.zst} .pvar).pvar.head.txt
tmp_pvar_body=${tmp_pvar_head%.head.txt}.body
tmp_pvar_split_prefix=$(basename ${tmp_pvar_body} | sed -e 's/original/split/g').

cat_or_zcat ${pvar} | egrep    '^#' > ${tmp_pvar_head}
cat_or_zcat ${pvar} | egrep -v '^#' > ${tmp_pvar_body}

split --lines 2000 --numeric-suffixes=1 --suffix-length 3 ${tmp_pvar_body} ${tmp_dir}/${tmp_pvar_split_prefix}

find ${tmp_dir} -type f -name "${tmp_pvar_split_prefix}*" | while read pvar_part ; do
    cat ${tmp_pvar_head} ${pvar_part} > ${data_d}/$(basename ${pvar_part}).pvar

    bash $(dirname ${SRCDIR})/helpers/pvar.to.vcf.sh ${data_d}/$(basename ${pvar_part}).pvar \
    | sed -e 's/chrMT/chrM/g' | sed -e 's/chrXY/chrX/g' > ${data_d}/$(basename ${pvar_part}).vcf
done

echo ${data_d}
