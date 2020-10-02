#!/bin/bash
set -beEuo pipefail

vep_d='/scratch/groups/mrivas/ukbb24983/cal/annotation_20201002'
n_batch=403

# out_f=${vep_d}/dev.ukb24983_cal_cALL_v2_hg19.vep101-loftee.tsv.gz
out_f=${vep_d}/ukb24983_cal_cALL_v2_hg19.vep101-loftee.tsv.gz

{
    cat ${vep_d}/output_vep/split.ukb24983_cal_cALL_v2_hg19.pvar.body.001.vep101-loftee.tsv | awk 'NR==1' | sed -e 's/^CHROM/#CHROM/g'

    seq -w ${n_batch} | while read idx_pad ; do
        f=${vep_d}/output_vep/split.ukb24983_cal_cALL_v2_hg19.pvar.body.${idx_pad}.vep101-loftee.tsv
        if [ -f ${f} ] ; then cat ${f} | awk 'NR>1' ; fi
    done

} | bgzip -l9 -@4 > ${out_f}

echo ${out_f}
