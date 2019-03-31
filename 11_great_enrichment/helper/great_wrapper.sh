#!/bin/bash
set -beEuo pipefail

usage () {
    echo "usage: $0 in.bed assembly out.dir"
}

get_onto_list () {
    assembly=$1
    onto_file="$(dirname $(dirname $(readlink -f $0)))/misc/great.ontology.${assembly}.20171029.txt"

    cat ${onto_file} | egrep -v '^#' | cut -f1 | tr "\n" "," | rev | cut -c2- | rev
}

if [ $# -lt 3 ] ; then usage >&2 ; exit 1 ; fi
in_file=$1
assembly=$2
out_d=$3

get_onto_list $assembly

exit 0
GREATER 

