#!/bin/bash
set -beEuo pipefail

cat_or_zcat () {
    local file=$1
    if [ ${file%.gz}.gz == ${file} ] ; then zcat ${file} ; else cat ${file} ; fi
}

phe_path_to_gwas_path () {
    local phe_path=$1
    local population=$2
    local gwas_type=$3
    echo ${phe_path} | sed -e "s%/phenotypedata/%/${gwas_type}/gwas/%g" | sed -e "s%/phe%/${population}%g"
}

GBE_ID_to_gwas_basename () {
    local GBE_ID=$1
    local gwas_variant_type=$2
    local prefix=$3
    echo ${prefix}.${GBE_ID}.${gwas_variant_type}.glm

#    ukb24983_v2_hg38.BIN1930.exome-spb.PHENO1.glm.logistic.hybrid.gz  # this Exome file name doesn't look right
#    ukb24983_v2_hg19.BIN10030500.genotyped.glm.logistic.hybrid.gz
}

get_dir_from_variant_type () {
    local variant_type=$1

    case ${variant_type} in
        "genotyped" ) 
            echo "cal" ;
            ;;
        "exome-spb.PHENO1" )
            echo "exome" ; # this is an adhock patch for Exome file names
            ;;
        "exome-spb" )
            echo "exome" ;
            ;;
        * )
            echo "" ;
            ;;
    esac
}

get_sumstats_name () {
    local phe_file=$1
    local population=$2
    local prefix=$3
    local variant_type=$4

    dir=$(phe_path_to_gwas_path $(dirname $phe_file) $population $(get_dir_from_variant_type ${variant_type}))
    basename_prefix=$(GBE_ID_to_gwas_basename $(basename ${phe_file%.phe}) ${variant_type} ${prefix})
    sumstats_file=$(find ${dir} -type f -name "${basename_prefix}*") 
    echo ${dir} >&2
    echo ${basename_prefix} >&2
    
    echo ${sumstats_file}
}

check_sumstats () {
    local phe_file=$1
    local population=$2
    local prefix=$3
    local variant_type=$4

    sumstats_file="$(get_sumstats_name ${phe_file} ${population} ${prefix} ${variant_type})"

    if [ ! -z ${sumstats_file} ] && [ -f ${sumstats_file} ] ; then
        wc_l=$(cat_or_zcat ${sumstats_file} | wc -l)
        wc_l_non_NAs=$(cat_or_zcat ${sumstats_file} | grep -v NA | wc -l)
    else
        sumstats_file="NA" ; wc_l="NA" ; wc_l_non_NAs="NA"
    fi
    echo "${sumstats_file} ${wc_l} ${wc_l_non_NAs}" | tr " " "\t"
}

check_sumstats_wrapper () {
    local phe_file_list=$1
    local idx=$2
    local population=$3
    local prefix=$4
    local variant_type=$5

    phe_file="$( cat ${phe_file_list} | awk -v idx=$idx 'NR==idx' )"
    echo "$idx $phe_file $(check_sumstats $phe_file $population $prefix ${variant_type})" | tr " " "\t"
}

