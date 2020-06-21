#!/bin/bash
set -beEuo pipefail

cd /oak/stanford/groups/mrivas/users/ytanigaw/sandbox/20200612_metal_power_analysis

for pop in NBW Afr SA EA rel others ; do
    bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_rg.sh --scratch LDSC_INI50_WB.sumstats.gz LDSC_INI50_${pop}.sumstats.gz LDSC_rg_WB_${pop}
done
