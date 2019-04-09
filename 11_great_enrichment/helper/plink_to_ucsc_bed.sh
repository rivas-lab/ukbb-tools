#!/bin/bash
set -beEuo pipefail

_default_n_top_hits=5000

usage () {
    echo "usage: $0 <in_file> [n_top_hits=${_default_n_top_hits}]"
}

cat_or_zcat () {
    file=$1
    if [ ${file%.gz}.gz == ${file} ] ; then zcat ${file} ; else cat ${file} ; fi
}

get_col_by_key () {
    file=$1
    key=$2
    cat_or_zcat $file | awk 'NR==1' | tr "\t" "\n" | awk -v key=$key '($0==key){print NR}'
}

get_thr () {
    file=$1
    key=$2
    num=$3
    cat_or_zcat $file | egrep -v '^#' | cut -f$( get_col_by_key $file $key ) \
    | grep -v "NA" | sort -n | awk -v num=$num 'NR==num'
}

plink_to_ucsc_bed () {
    file=$1
    key=$2
    num=$3
    cat_or_zcat $file | egrep -v '^#' \
    | cut -f1-3,$( get_col_by_key $file $key ) \
    | sed -e 's/^23/X/g' | sed -e 's/^24/Y/g' | sed -e 's/^25/X/g' | sed -e 's/^26/MT/g' | awk '{print "chr" $0}'\
    | awk -v thr=$( get_thr $file $key $num ) -v na="NA" -v OFS='\t' '($4 <= thr && $4 != na){print $1, $2, $2 + 1, $3}'
}

if [ $# -lt 1 ] ; then usage >&2 ; exit 1; fi 
in_file=$1
if [ $# -gt 1 ] ; then n_top_hits=$2 ; else n_top_hits=${_default_n_top_hits} ; fi

# sort by P-value
plink_to_ucsc_bed $in_file "P" $n_top_hits

