#!/bin/bash
# set -beEuo pipefail

show_master_file_header () {
cat <<- EOF
	# === DESCRIBE AND PROCESS THE INPUT INPUT FILES ===
	
	SCHEME STDERR
	MARKERLABEL ID
	ALLELE A1 REF
	EFFECT BETA
	PVALUE P
	STDERR SE
	WEIGHT OBS_CT
	
	EFFECT_PRINT_PRECISION 7
	STDERR_PRINT_PRECISION 7
	
EOF
}

show_master_file_body () {
    local file_list=$1
    cat ${file_list} \
    | awk -v prefix="PROCESS " '(length($0)>0){print prefix $0}'
}

show_master_file_footer () {
    local out_file=$1
cat <<- EOF
	
	OUTFILE ${out_file} .tbl
	MINWEIGHT 10000
	
	ANALYZE HETEROGENEITY
	
	QUIT
EOF
}

show_master_file () {
    local file_list=$1
    local out_file=$2

    show_master_file_header

    show_master_file_body ${file_list}

    show_master_file_footer ${out_file}
}

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

show_header () {
    local file=$1
    # ! cat_or_zcat $file | head -n1
    cat_or_zcat $file | egrep '^#'
}

get_col_idx () {
    # Return the column index
    local file=$1
    local key=$2
    show_header $file | sed -e "s/^#//g" | tr "\t" "\n" | awk -v key=$key '($0 == key){print NR}'
}

add_BETA_from_OR () {
    # For a given summary statistics (with OR column), we compute BETA as log(OR)
    # We also replace LOG(OR)_SE with SE in the header line
    local in_file=$1

    local col_OR=$( get_col_idx $in_file "OR")

    echo "$(show_header $in_file) BETA" | sed -e "s/LOG(OR)_SE/SE/g" | tr " " "\t"

    cat_or_zcat ${in_file} \
    | egrep -v '^#' \
    | awk -v OFS='\t' -v cOR=${col_OR} '(NR==1 || $cOR != "NA"){print $0, log($cOR)}'
}

metal_pre_processing () {
    # We apply pre-processing
    #   1. For sumstats from logistic regression (we detect OR column), we add BETA column
    #   2. We focus on the lines where P != "NA" and ERRCODE == "." (if ERRCODE column exists)
    local in_file=$1
    local ERRCODE_filter="FALSE"
    if [ $# -gt 1 ] && [ $2 == "ERRCODE" ] ; then ERRCODE_filter="TRUE" ; fi

    check_OR_flag=$(show_header $in_file | tr "\t" "\n" | cat /dev/stdin <(echo OR) | grep OR | wc -l)
    check_ERRCODE_flag=$(show_header $in_file | tr "\t" "\n" | cat /dev/stdin <(echo ERRCODE) | grep ERRCODE | wc -l)
    if [ "${ERRCODE_filter}" == "TRUE" ] && [ "${check_ERRCODE_flag}" -eq 1 ] ; then
        echo "ERRCODE filter option is enabled but the ERRCODE column does not exist in the input file: ${in_file}" >&2
        exit 1
    fi   
    
    col_P=$( get_col_idx $in_file "P")
    if [ "${ERRCODE_filter}" == "TRUE" ] && [ "${check_ERRCODE_flag}" -gt 1 ] ; then
        col_ERRCODE=$( get_col_idx $in_file "ERRCODE")
    fi

    if [ "${check_OR_flag}" -gt 1 ] ; then
        add_BETA_from_OR ${in_file}
    else
        cat_or_zcat ${in_file}
    fi | if [ "${ERRCODE_filter}" == "TRUE" ] && [ "${check_ERRCODE_flag}" -gt 1 ] ; then
        awk -v OFS='\t' -v cP=${col_P} -v cE=${col_ERRCODE} '((NR == 1) || ($cP != "NA" && $cE == "."))'
    else
        awk -v OFS='\t' -v cP=${col_P} '((NR == 1) || ($cP != "NA"))'
    fi
}

show_metal_input_files () {
    local metal_info=$1
    nr_s=$(cat_or_zcat ${metal_info} | egrep -n '^# == original input files ==' | awk -v FS=':' '{print $1}')
    nr_e=$(cat_or_zcat ${metal_info} | egrep -n '^# == METAL info file ==' | awk -v FS=':' '{print $1}')

    cat_or_zcat ${metal_info} | awk -v nr_s=${nr_s}  -v nr_e=${nr_e} 'nr_s < NR && NR < nr_e'
}

