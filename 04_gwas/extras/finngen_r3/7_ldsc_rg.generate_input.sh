#!/bin/bash
set -beEuo pipefail

finngen_d="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3"
ukb_WB_d="/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/white_british"
finngen_idx="/oak/stanford/groups/mrivas/public_data/finngen_r3/ldsc_h2.tsv"
ukb_idx="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/05_gbe/phenotype_info.tsv"
ukb_h2="/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/h2.white_british.tsv"
ukb_WB_min_N=1000

out_prefix=$1
out_ukb="${out_prefix}.ukb.tsv"
out_finngen="${out_prefix}.finngen.tsv"

{
    echo "#GBE_ID WB_f"

    cat ${ukb_idx} | awk -v ukb_WB_min_N=${ukb_WB_min_N} -v FS='\t' '(NR>1 && $8 >= ukb_WB_min_N){print $1}' | egrep '^HC' | sort | comm -12 /dev/stdin <(cat ${ukb_h2} | egrep '^HC' |  awk '($2 > 0){print $1}' | sort) | sort | comm -12 /dev/stdin <(find ${ukb_WB_d} -type f -name "*.sumstats.gz" | awk -v FS='/' '{print $NF}' | awk -v FS='.' '{print $2}' | sort) | while read GBE_ID ; do
        UKB_f="${ukb_WB_d}/ukb24983_v2_hg19.${GBE_ID}.array-combined.sumstats.gz"
        echo "${GBE_ID} ${UKB_f}"
    done
}  | tr ' ' '\t' > ${out_ukb}

{
    echo "#FinnGen FinnGen_f"

    cat ${finngen_idx} | awk '(NR>1 && $2 > 0){print $1}' | while read FinnGen ; do

        FinnGen_f="${finngen_d}/summary_stats_hg19_ldsc_munge/finngen_r3_${FinnGen}.hg19.sumstats.gz"

        echo "${FinnGen} ${FinnGen_f}"
    done
} | tr ' ' '\t' > ${out_finngen}

echo ${out_prefix}
echo ${out_ukb}
echo ${out_finngen}
