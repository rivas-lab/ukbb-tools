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
    #   2. We remove the lines where P = NA
    local in_file=$1

    col_P=$( get_col_idx $in_file "P")
    check_OR_flag=$(show_header $in_file | tr "\t" "\n" | cat /dev/stdin <(echo OR) | grep OR | wc -l)

    if [ ${check_OR_flag} -gt 1 ] ; then
        echo "applying add_BETA_from_OR to ${in_file}"
        add_BETA_from_OR ${in_file}
    else
        cat_or_zcat ${in_file}
    fi | awk -v OFS='\t' -v cP=${col_P} '$cP != "NA"'
}

extract_loci () {
    # For a given input file, we extract the list of loci (CHROM, POS, and ID)
    local in_file=$1
    
    local col_CHROM=$( get_col_idx $in_file "CHROM")
    local col_POS=$(   get_col_idx $in_file "POS")
    local col_ID=$(    get_col_idx $in_file "ID")
   
    cat_or_zcat ${in_file} \
    | awk -v OFS='\t' -v cCHROM=${col_CHROM} -v cPOS=${col_POS} -v cID=${col_ID} \
    '(NR>1){print $cCHROM, $cPOS, $cID}'
}

extract_loci_for_files () {
    # For a given set of files, we extract the loci (CHROM, POS, and ID)
    local in_files=$1 # a file containing the list of file
    if [ $# -gt 1 ] ; then nCores=$2 ; else nCores=1 ; fi
    
    echo "#CHROM POS ID" | tr " " "\t"
    
    cat_or_zcat ${in_files} | while read f ; do 
        extract_loci $f
    done | sort --parallel ${nCores} -k1,1V -k2,2n -k3,3 -u
}
