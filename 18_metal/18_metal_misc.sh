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
    local outfile=$1
cat <<- EOF
	
	OUTFILE ${outfile} .tbl
	MINWEIGHT 10000
	
	ANALYZE HETEROGENEITY
	
	QUIT
EOF
}

show_master_file () {
    local file_list=$1
    local outfile=$2

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
    local file=$1
    local key=$2
    show_header $file | sed -e "s/^#//g" | tr "\t" "\n" | awk -v key=$key '($0 == key){print NR}'
}

add_BETA_from_OR () {
    local in_file=$1

    local col_OR=$( get_col_idx $in_file "OR")

    echo "$(show_header $in_file) BETA" | tr " " "\t"

    cat_or_zcat ${in_file} \
    | egrep -v '^#' \
    | awk -v OFS='\t' -v cOR=${col_OR} '{print $0, log($cOR)}'
}

extract_loci () {
    local in_file=$1
    
    local col_CHROM=$( get_col_idx $in_file "CHROM")
    local col_POS=$(   get_col_idx $in_file "POS")
    local col_ID=$(    get_col_idx $in_file "ID")
   
    cat_or_zcat ${in_file} \
    | awk -v OFS='\t' -v cCHROM=${col_CHROM} -v cPOS=${col_POS} -v cID=${col_ID} \
    '(NR>1){print $cCHROM, $cPOS, $cID}'
}

extract_loci_for_files () {
    local in_files=$1
    if [ $# -gt 1 ] ; then nCores=$2 ; else nCores=1 ; fi
    
    echo "#CHROM POS ID" | tr " " "\t"
    
    cat_or_zcat ${in_files} | while read f ; do 
        extract_loci $f
    done | sort --parallel ${nCores} -k1,1V -k2,2n -k3,3 -u
}