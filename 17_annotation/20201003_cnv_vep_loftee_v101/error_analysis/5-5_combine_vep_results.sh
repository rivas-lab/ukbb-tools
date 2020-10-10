#!/bin/bash
set -beEuo pipefail

vep_d='/scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003'

unfinished_v=20201007
n_batch=5040

out_f=${vep_d}/dev.cnv.${unfinished_v}.unfinished.vep101-loftee.tsv.gz

{
    cat ${vep_d}/output_vep_${unfinished_v}/split.cnv.${unfinished_v}.unfinished.pvar.body.${n_batch}.vep101-loftee.tsv | awk 'NR==1' | sed -e 's/^CHROM/#CHROM/g'

    seq -w ${n_batch} | while read idx_pad ; do
        f=${vep_d}/output_vep_${unfinished_v}/split.cnv.${unfinished_v}.unfinished.pvar.body.${idx_pad}.vep101-loftee.tsv
        if [ -f ${f} ] ; then cat ${f} | awk 'NR>1' ; fi
    done

} | bgzip -l9 -@4 > ${out_f}

echo ${out_f}

