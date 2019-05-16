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
    echo ${prefix}.${GBE_ID}.${gwas_variant_type}.PHENO1.glm

#    ukb24983_v2_hg38.BIN1930.exome-spb.PHENO1.glm.logistic.hybrid.gz  
#    ukb24983_v2_hg19.BIN10030500.genotyped.PHENO1.glm.logistic.hybrid.gz
}

get_dir_from_variant_type () {
    local variant_type=$1

    case ${variant_type} in
        "genotyped" ) 
            echo "cal" ;
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

create_colnames_tmpfile () {
    local file=$1
    if [ $# -gt 1 ] ; then tmp_dir_root=$2 ; else tmp_dir_root=/dev/shm ; fi
    tmp_colnames=$(mktemp -p $tmp_dir_root); 
    cat_or_zcat $file | awk 'NR==1' | sed -e 's/^#//g' | tr "\t" "\n" > ${tmp_colnames}
    echo ${tmp_colnames}
}

find_col_idx_from_colnames_file () {
    local colnames_file=$1
    local key=$2

#    cat $colnames_file | awk -v key=$key '(key==key){print NR}'
    cat $colnames_file | egrep -n $key | awk -v FS=':' '{print $1}' | tr "\n" "," | rev | cut -c2- | rev
}


find_col_idxs_from_colnames_file () {
    local colnames_file=$1
    shift 1
    for key in $@ ; do
        find_col_idx_from_colnames_file $colnames_file $key
    done | tr "\n" "," | sed -e 's/,$//g'
}

show_sumstats_header () {
    echo "CHROM POS Variant_ID GBE_ID REF ALT A1 BETA SE P STAT" | tr " " "\t" 
}


show_sumstats () {
    local GBE_ID=$1
    local sumstats_file=$2
    if [ $# -gt 2 ] ; then p_val_thr=$3 ; else p_val_thr=1 ; fi

    if [ -f $sumstats_file ] ; then
        local cols_tmp=$(create_colnames_tmpfile $sumstats_file)
        local cols=$(find_col_idxs_from_colnames_file $cols_tmp 'CHROM' 'POS' 'ID' 'REF' 'ALT' 'A1' 'BETA|OR' 'SE' 'Z_STAT|T_STAT' '^P$')
        rm $cols_tmp
        cat_or_zcat $sumstats_file | awk 'NR>1' | grep -v 'NA' | cut -f $cols \
            | awk -v OFS='\t' -v GBE_ID=${GBE_ID} -v p_val_thr=${p_val_thr} '($10 <= p_val_thr){print $1, $2, $3, GBE_ID, $4, $5, $6, $7, $8, $10, $9}'
    fi
}

