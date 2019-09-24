#!/bin/bash
set -beEuo pipefail

pops=(
"white_british"
"non_british_white"
"african"
"s_asian"
"e_asian"
)

check_ldmap_out () {
    local default_out_d="/scratch/groups/mrivas/ukbb/24983/imp/ldmap"

    local c="$1"
    local pop="$2"
    if [ $# -gt 2 ] ; then out_d=$3 ; else out_d=${default_out_d} ; fi

    for ext in "bool.prune.in" "bool.orune.out" "ld_map.ld.gz" ; do
        local file="${out_d}/ukb24983_imp_chr${c}_v3.${pop}.${ext}"
        if [ ! -f ${file} ] ; then
            echo "${c} ${pop} ${ext} ${file}" | tr " " "\t"
        fi
    done
}

echo "The following files are missing.." >&2

for c in $(seq 1 22) X XY ; do for pop in "${pops[@]}" ; do
    check_ldmap_out $c $pop
done; done

