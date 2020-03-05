#!/bin/bash
set -beEuo pipefail

pops=(
"white_british"
"non_british_white"
"african"
"s_asian"
"e_asian"
)

source /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/14_LD_map/14_LD_map_misc.sh

for pop in ${pops[@]} ; do
    out_prefix="/scratch/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/ldmap/ukb24983_cal_hla_imp.${pop}"
    check_ldmap_out ${out_prefix}
done
