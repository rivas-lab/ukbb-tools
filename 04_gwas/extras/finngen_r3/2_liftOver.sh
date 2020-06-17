#!/bin/bash
set -beEuo pipefail

data_d="/scratch/groups/mrivas/public_data/summary_stats/finngen_r3"

# input
sumstats_with_largest_N="${data_d}/summary_stats/finngen_r3_M13_MUSCULOSKELETAL.gz"

# output
out_hg38=${data_d}/finngen_r3_variants.tsv.gz
out_hg19=${data_d}/finngen_r3_variants.hg19
out_hg19_unmapped=${data_d}/finngen_r3_variants.hg19.unmapped
out_hg19_flipcheck=${data_d}/finngen_r3_variants.hg19.fasta.tsv.gz

echo "extracting the list of variants..."

{
    echo "#CHROM POS REF ALT ID rsids nearest_genes" | tr " " "\t"
    zcat ${sumstats_with_largest_N} \
    | awk -v sep=':' -v FS='\t' -v OFS='\t' '(NR>1){ print $1, $2, $3, $4, $1 sep $2 sep $3 sep $4, $5, $6 }' 
} | bgzip -l9 -@6 > ${out_hg38}

echo ${out_hg38}

echo "running liftOver ..."

bash ~/repos/rivas-lab/ukbb-tools/09_liftOver/liftOver_wrapper.sh --src_genome hg38 --dst_genome hg19 ${out_hg38} ${out_hg19} ${out_hg19_unmapped}

echo "fetching the ref allele from FASTA file ..."

bash ~/repos/rivas-lab/ukbb-tools/09_liftOver/flipcheck.sh ${out_hg19}.gz | bgzip -l9 -@6 > ${out_hg19_flipcheck}

rm ${out_hg19}.gz

echo ${out_hg19_unmapped}
echo ${out_hg19_flipcheck}
