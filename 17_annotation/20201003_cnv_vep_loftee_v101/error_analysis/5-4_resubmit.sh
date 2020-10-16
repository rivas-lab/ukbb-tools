#!/bin/bash
set -beEuo pipefail

idx=$1
list=$2
if [ $# -gt 2 ] ; then offset=$3 ; else offset=0 ; fi

idxx=$(cat $list | egrep -v '^#' | awk -v l=$(perl -e "print($idx + $offset)") 'NR==l' | sed -e 's/^0\+//g')

bash 5-2_run_vep.sh $idxx
