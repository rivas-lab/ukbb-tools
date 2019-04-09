#!/bin/bash
set -beEuo pipefail

usage () {
    echo "usage: $0 in.bed assembly out.dir"
}

source $(dirname $(readlink -f $0))/great_misc_func.sh

if [ $# -lt 3 ] ; then usage >&2 ; exit 1 ; fi
in_file=$1
assembly=$2
out_d=$3

GREATER \
    --requiredTests=neither \
    --ontologies=$( get_onto_list $assembly ) \
    $assembly $in_file $out_d

