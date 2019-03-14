#!/bin/bash
set -beEuo pipefail

usage () {
    echo "extract phenotype" >&2
    echo "usage: ${0} GBE_ID [covar] [master_phe_file]" >&2
}

default_master_phe="$OAK/dev-ukbb-tools/phewas/resources/master.phe"

if [ $# -lt 1 ] ; then usage ; exit 1 ; fi
GBE_ID=$1
if [ $# -gt 2 ] ; then master_phe=$3 ; else master_phe=${default_master_phe} ; fi 

phe_col_id=$( cat $master_phe | awk 'NR==1' | tr "\t" "\n" | egrep -n "${GBE_ID}$" | awk -v FS=':' '{print $1}' )

if [ $# -gt 1 ] && [ $2 == 'covar' ] ; then
    cat $master_phe | cut -f1-45,${phe_col_id}
else
    cat $master_phe | cut -f1-2,${phe_col_id}
fi

