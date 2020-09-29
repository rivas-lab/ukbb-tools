#!/bin/bash
set -beEuo pipefail

cat_or_zstdcat () {
    local file=$1
    if [ "${file%.zst}.zst" == "${file}" ] ; then zstdcat $file ; else cat $file ; fi
}

pgen_to_bed () {
    local pvarfile=$1
    local bfile=$2
    local plink_mem=$3
    local cores=$4
    local plink_cmd=$5
    
    local pfile=$(echo ${zstfile%.zst} | sed -e "s/.pvar$//g")    
    
    echo plink2 --memory ${plink_mem} --threads ${cores} --pfile "${pfile}" $([ "${pvarfile%.zst}.zst" == "${pvarfile}" ] && echo "vzs" || echo "") --out "${bfile}" --${plink_cmd}
    plink2      --memory ${plink_mem} --threads ${cores} --pfile "${pfile}" $([ "${pvarfile%.zst}.zst" == "${pvarfile}" ] && echo "vzs" || echo "") --out "${bfile}" --${plink_cmd}
}

prep_files () {
    local zstfile=$1
    local prefix=$2
    local tmp_dir=$3
    local helper_R_script=$4
    local plink_mem=$5
    local cores=$6
    
    local pfile=$(echo $zstfile | sed -e "s/.zst$//g" | sed -e "s/.pvar$//g")
    local tmp_pvar="${tmp_dir}/${prefix}.tmp.in.pvar"
    local tmp_bim="${tmp_dir}/$(basename $pfile).bim"
    local tmp_bed="${tmp_dir}/$(basename $pfile).bed"
    local tmp_fam="${tmp_dir}/$(basename $pfile).fam"

    echo ${tmp_bed}
    echo ${tmp_bim}
    echo ${tmp_fam}

    # generate bed/fam file from pfile when needed.
    # we will generate bim file later
    if [ ! -s ${tmp_bed} ] ; then
        # bed file does not exist --> generate bed + fam (+ bim) file
        pgen_to_bed ${zstfile} ${tmp_bed%.bed} ${plink_mem} ${cores} "make-bed"
        mv ${tmp_bed%.bed}.log ${tmp_bed%.bed}.bed.log
    elif [ ! -s ${tmp_fam} ] ; then
        # fam file does not exist --> generate fam (+ bim) file
        pgen_to_bed ${zstfile} ${tmp_bed%.bed} ${plink_mem} ${cores} "make-just-fam"
        mv ${tmp_bed%.bed}.log ${tmp_bed%.bed}.fam.log
    fi
    
    # generate the short variant ID file
    if [ ! -s ${pfile}.IDs.tsv ] ; then
        cat_or_zstdcat ${zstfile} | grep -v '##' > ${tmp_pvar}
        Rscript ${helper_R_script} ${tmp_pvar} ${prefix} ${pfile}.IDs.tsv
        rm ${tmp_pvar}
    fi
    
    # write to a bim file
    cat ${pfile}.IDs.tsv \
    | awk -v FS='\t' -v OFS='\t' '(NR>1){print $1, $6, 0, $2, $5, $4}' > ${tmp_bim}
}

