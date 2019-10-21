#!/bin/bash
set -beEuo pipefail

chr=$1
if [ $# -gt 1 ] ; then cpu=$2 ; else cpu=10 ; fi
if [ $# -gt 2 ] ; then mem=$3 ; else mem=150000 ; fi

tmp_dir_root=${LOCAL_SCRATCH:=/tmp}/u/$USER
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

assign_var_id () {
    local pvar_file=$1
    local tmp_dir=$2   

    local tmp_pvar=${tmp_dir}/$(basename ${pvar_file%.zst})
    local info_id="ORIGINAL_VAR_ID"

    #zstdgrep '#' ${pvar_file} >> ${tmp_pvar}
    echo "##INFO=<ID=${info_id},Number=1,Type=String,Description=\"The variant ID in the original data release\">" > ${tmp_pvar}
    echo "#CHROM POS ID REF ALT INFO" | tr " " "\t" >> ${tmp_pvar}

    zstdgrep -v '#' ${pvar_file} \
        | awk -v OFS='\t' -v sep=":" -v info_id=${info_id} \
        '{print $1, $2, $1 sep $2 sep $4 sep $5, $4, $5, info_id "=" $3}' \
        | sed -e "s/.;${info_id}=.$/./g" \
        | sed -e "s/.;${info_id}=/${info_id}=/g" >> ${tmp_pvar}
    
    zstd --rm -9 ${tmp_pvar}
    mv ${tmp_pvar}.zst ${pvar_file}
}

pvar_file=/oak/stanford/groups/mrivas/ukbb/24983/imp/pgen/ukb24983_imp_chr${chr}_v3.pvar.zst
cp $pvar_file ${pvar_file%.pvar.zst}.original.pvar.zst
cp ${pvar_file%.pvar.zst}.bim ${pvar_file%.pvar.zst}.original.bim
assign_var_id ${pvar_file} ${tmp_dir}
plink2 --threads ${cpu} --memory ${mem} --pfile ${pvar_file%.pvar.zst} vzs --make-just-bim --out ${pvar_file%.pvar.zst}

