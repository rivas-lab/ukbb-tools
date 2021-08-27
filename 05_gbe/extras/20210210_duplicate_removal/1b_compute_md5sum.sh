#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

source 0_parameters.sh


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


get_col_idx () {
    local phe=$1
    local phe_ID=$2
    
    ! cat_or_zcat ${phe} | head -n1 | tr '\t' '\n' \
    | awk -v ID=${phe_ID} '($0==ID){print NR}'
} 

compute_md5sum () {
    local master_phe_f=$1
    local GBE_ID=$2
    
    idx=$(get_col_idx ${master_phe_f} ${GBE_ID})
    cat_or_zcat ${master_phe_f} | cut -f${idx} | tail -n+2 | md5sum \
    | awk -v ID=${GBE_ID} -v OFS='\t' '{print ID, $1}'

}

# copy the phe file to SSD
ls_m_phe=${tmp_dir}/$(basename ${tmp_master_phe})
cp ${tmp_master_phe} ${ls_m_phe}

{
    echo "#GBE_ID md5sum" | tr ' ' '\t'
    cat ${potential_dups_GBE_f} | awk '(NR>1){print $1}' | while read GBE_ID ; do

        compute_md5sum ${ls_m_phe} ${GBE_ID}

    done
} > ${potential_dups_GBE_md5_f}
