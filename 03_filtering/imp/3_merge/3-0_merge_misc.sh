#!/bin/bash
set -beEuo pipefail

cat_or_zstdcat () {
    local file=$1
    if [ "${file%.zst}.zst" == "${file}" ] ; then zstdcat $file ; else cat $file ; fi
}

prep_files () {
    local zstfile=$1
    local prefix=$2
    local tmp_dir=$3
    local helper_R_script=$4
    
    local pfile=$(echo $zstfile | sed -e "s/.zst$//g" | sed -e "s/.pvar$//g")
    local tmp_pvar="${tmp_dir}/${prefix}.tmp.in.pvar"
    local tmp_bim="${tmp_dir}/$(basename $pfile).bim"
    local tmp_bed="${tmp_dir}/$(basename $pfile).bed"
    local tmp_fam="${tmp_dir}/$(basename $pfile).fam"

    ln -s ${pfile}.fam ${tmp_fam}    
    ln -s ${pfile}.bed ${tmp_bed}
    
    cat_or_zstdcat ${zstfile} | grep -v '##' > ${tmp_pvar}
    Rscript ${helper_R_script} ${tmp_pvar} ${prefix} ${pfile}.IDs.tsv
    rm ${tmp_pvar}
    
    cat ${pfile}.IDs.tsv \
    | awk -v FS='\t' -v OFS='\t' '(NR>1){print $1, $6, 0, $2, $5, $4}' > ${tmp_bim}
}
