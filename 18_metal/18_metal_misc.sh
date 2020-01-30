#!/bin/bash
set -beEuo pipefail

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
