#!/bin/bash
set -beEuo pipefail

_default_n_top_hits=5000

usage () {
    echo "usage: $0 <in_file> [n_top_hits=${_default_n_top_hits}]"
}

source $(dirname $(readlink -f $0))/great_misc_func.sh

get_thr_clumped () {
    local file=$1
    local num=$2
    cat_or_zcat $file | awk '(NR>1 && length($0)> 0){print $5}'  \
    | sort -g | awk -v num=$num 'NR==num'
}

clumped_to_ucsc_bed () {
    local file=$1
    local num=$2

    local wc_l=$( cat_or_zcat $file | awk '(NR>1 && length($0) > 0)' | wc -l )
    local thr=$(get_thr_clumped $file $num)
    cat_or_zcat $file \
        | awk '(NR>1 && length($0) > 0){print $1, $4, $3, $5}' \
    | sed -e 's/^23/X/g' | sed -e 's/^24/Y/g' | sed -e 's/^25/X/g' | sed -e 's/^26/MT/g' | awk '{print "chr" $0}'\
    | if [ $num -lt $wc_l ] ; then
        awk -v thr=$thr -v na="NA" -v OFS='\t' '($4 <= thr && $4 != na){print $1, $2, $2 + 1, $3 "_" $4}'
    else
        awk -v na="NA" -v OFS='\t' '($4 != na){print $1, $2, $2 + 1, $3 "_" $4}'
    fi
}

if [ $# -lt 1 ] ; then usage >&2 ; exit 1; fi 
in_file=$1

if [ $# -gt 1 ] ; then n_top_hits=$2 ; else n_top_hits=${_default_n_top_hits} ; fi

# sort by P-value
clumped_to_ucsc_bed $in_file $n_top_hits

