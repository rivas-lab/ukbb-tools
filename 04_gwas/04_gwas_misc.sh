#!/bin/bash
#set -beEuo pipefail

#########################################################
# Common helper functions
#########################################################

cat_or_zcat () {
    local file=$1
    if [ ${file%.gz}.gz == ${file} ] ; then zcat ${file} ; else cat ${file} ; fi
}

#########################################################
# Helper functions used in gwas.py
#########################################################

find_file () {
    local filePrefix=$1
    find $(dirname $filePrefix) -name "$(basename $filePrefix)*glm*"
}

find_suffix () {
    local filePrefix=$1
    find_file $filePrefix | sed -e "s%${filePrefix}.%%"
}

mv_glm_file_and_bgzip () {
    echo "Moving glm file and bgzipping..."
    local filePrefix=$1
    local src=$(find_file $filePrefix | egrep -v 'gz$')
    local dst=$filePrefix.$(find_suffix $filePrefix | egrep -v 'gz$' | cut -d'.' -f2-)
    echo "Source: ${src}"
    echo "Destination: ${dst}"
    mv $src $dst
    apply_bgzip $dst
}

apply_bgzip () {
    local file=$1
    if [ "${file%.gz}.gz" != "${file}" ] ; then bgzip -l9 -f ${file} ; fi
}

get_log_filename () {
    local filePrefix=$1
    
    echo "$(dirname ${filePrefix})/logs/$(basename ${filePrefix}).log"
}

post_processing () {
    local filePrefix=$1
    local logFile=$(get_log_filename ${filePrefix})
    if [ ! -d $(dirname ${logFile}) ] ; then mkdir -p $(dirname ${logFile}) ; fi
    echo "Moving the log file..."
    echo "Source: ${filePrefix}.log"
    echo "Destination: ${logFile}"
    mv ${filePrefix}.log ${logFile}
    mv_glm_file_and_bgzip $filePrefix 
}

combine_two_sumstats () {
# combine two summary statistics (as well as the corresponding log files)
    local inFile1Prefix=$1
    local inFile2Prefix=$2
    local outFilePrefix=$3
    echo "inFile1Prefix = $inFile1Prefix"
    echo "infile2Prefix = $inFile2Prefix"
    echo "outFilePrefix = $outFilePrefix"
    
    if [ $# -gt 3 ] ; then threads=$4 ; else threads=1 ; fi
    
    local file1=$(find_file $inFile1Prefix | egrep 'gz$')
    local file2=$(find_file $inFile2Prefix | egrep 'gz$')
    local ending=$(find_suffix $inFile1Prefix | egrep 'gz$')

    echo "file1 = $file1"
    echo "file2 = $file2"
    echo "ending = $ending"
    cat $(get_log_filename ${inFile1Prefix}) $(get_log_filename ${inFile2Prefix}) > $(get_log_filename ${outFilePrefix})
    rm $(get_log_filename ${inFile1Prefix}) $(get_log_filename ${inFile2Prefix})

    { 
    zcat $file1 | egrep '^#' #header
    zcat $file1 $file2 | egrep -v '^#' | sort -k1,1V -k2,2n -u --parallel=${threads}
    } | bgzip -l9 -f > ${outFilePrefix}.${ending%.gz}.gz
    rm $file1 $file2
}

#########################################################
# Common helper functions
#########################################################

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
    local pop=$4
    if [ $pop = "white_british" ]; then
        echo ${prefix}.${GBE_ID}.${gwas_variant_type}.PHENO1.glm
    else
        echo ${prefix}.${GBE_ID}.${gwas_variant_type}.glm
    fi
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
    basename_prefix=$(GBE_ID_to_gwas_basename $(basename ${phe_file%.phe}) ${variant_type} ${prefix} ${population})
    sumstats_file=$(find ${dir} -type f -name "${basename_prefix}*gz") 
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

        local col_beta="$(find_col_idx_from_colnames_file $cols_tmp 'BETA')"
        local col_or="$(find_col_idx_from_colnames_file $cols_tmp 'OR')"

        rm $cols_tmp

        if [ "${col_beta}" != "" ] && [ "${col_or}" == "" ] ; then
            cat_or_zcat $sumstats_file | awk 'NR>1' | grep -v 'NA' | cut -f $cols \
                | awk -v OFS='\t' -v GBE_ID=${GBE_ID} -v p_val_thr=${p_val_thr} '($10 <= p_val_thr){print $1, $2, $3, GBE_ID, $4, $5, $6, $7, $8, $10, $9}' 
                
        elif [ "${col_or}" != "" ] && [ "${col_beta}" == "" ] ; then
            cat_or_zcat $sumstats_file | awk 'NR>1' | grep -v 'NA' | cut -f $cols \
                | awk -v OFS='\t' -v GBE_ID=${GBE_ID} -v p_val_thr=${p_val_thr} '($10 <= p_val_thr){print $1, $2, $3, GBE_ID, $4, $5, $6, log($7), $8, $10, $9}'
        fi
    fi       
}


get_field_from_pop () {
    local pop=$1
    local const=7
    echo "white_british non_british_white african e_asian s_asian" \
        | tr " " "\n" | awk -v pop=$pop -v const=$const '($0 == pop){print NR + const}'
}

