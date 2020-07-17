#!/bin/bash
set -beEuo pipefail

wc_l () {
    local f=$1

    echo $(zcat $f | wc -l ) $f | tr ' ' '\t'
}

export -f wc_l

{
echo "#wc_l f" | tr ' ' '\t'
find data/ -name "*.gz" | sort | parallel -k wc_l
} | tee 4_wc_l.$(date +%Y%m%d-%H%M%S).tsv
