#!/bin/bash
set -beEuo pipefail

cd /oak/stanford/groups/mrivas/users/ytanigaw/sandbox/20200612_metal_power_analysis

bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/white_british/ukb24983_v2_hg19.INI50.array-combined.glm.linear.gz LDSC_INI50_WB

bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/non_british_white/ukb24983_v2_hg19.INI50.array-combined.glm.linear.gz LDSC_INI50_NBW

bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/african/ukb24983_v2_hg19.INI50.array-combined.glm.linear.gz LDSC_INI50_Afr

bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/s_asian/ukb24983_v2_hg19.INI50.array-combined.glm.linear.gz LDSC_INI50_SA

bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/e_asian/ukb24983_v2_hg19.INI50.array-combined.glm.linear.gz LDSC_INI50_EA

bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch /oak/stanford/groups/mrivas/projects/related/ukb24983_v2_hg19.1915,1925,1940,1944,1954,1962,1964,1965,1966,1967.array-combined.INI50.glm.linear.gz LDSC_INI50_rel

bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch /oak/stanford/groups/mrivas/projects/gwas_others/ukb24983_v2_hg19.1915,1925,1940,1944,1954,1962,1964,1965,1966,1967.array-combined.INI50.glm.linear.gz LDSC_INI50_others
