#!/bin/bash
set -beEuo pipefail

vep_d='/scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003'
n_batch=998

out_f=${vep_d}/dev.cnv.vep101-loftee.tsv.gz
# out_f=${vep_d}/cnv.vep101-loftee.tsv.gz

{
    cat ${vep_d}/output_vep/split.cnv.pvar.body.001.vep101-loftee.tsv | awk 'NR==1' | sed -e 's/^CHROM/#CHROM/g'

    seq -w ${n_batch} | while read idx_pad ; do
        f=${vep_d}/output_vep/split.cnv.pvar.body.${idx_pad}.vep101-loftee.tsv
        if [ -f ${f} ] ; then cat ${f} | awk 'NR>1' ; fi
    done

} | bgzip -l9 -@4 > ${out_f}

echo ${out_f}
