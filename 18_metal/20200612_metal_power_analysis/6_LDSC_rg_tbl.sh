#!/bin/bash
set -beEuo pipefail

cd /oak/stanford/groups/mrivas/users/ytanigaw/sandbox/20200612_metal_power_analysis

{
bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_rg_view.sh LDSC_rg_WB_NBW.log
for pop in Afr SA EA rel others ; do
    bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_rg_view.sh LDSC_rg_WB_${pop}.log | awk 'NR>1'
done
} | sed -e "s%/oak/stanford/groups/mrivas/users/ytanigaw/sandbox/20200612_metal_power_analysis/%%g" | sed -e "s/.sumstats.gz//g" | sed -e "s/LDSC_//g" > $(dirname $(readlink -f $0))/6_LDSC_rg_tbl.tsv
