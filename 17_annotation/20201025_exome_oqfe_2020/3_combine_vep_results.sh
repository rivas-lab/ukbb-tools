#!/bin/bash
set -beEuo pipefail

vep_d='/scratch/groups/mrivas/ukbb24983/exome-oqfe2020-annotation'
n_batch=900

out_f=${vep_d}/UKBexomeOQFE.vep101.tsv.gz

{
    cat ${vep_d}/output_vep/split.UKBexomeOQFE.pvar.body.001.vep101-loftee.tsv | awk 'NR==1' | sed -e 's/^CHROM/#CHROM/g'

    seq -w ${n_batch} | while read idx_pad ; do
        f=${vep_d}/output_vep/split.UKBexomeOQFE.pvar.body.${idx_pad}.vep101-loftee.tsv 
        if [ -f ${f} ] ; then cat ${f} | awk 'NR>1' ; fi
    done

} | bgzip -l9 -@4 > ${out_f}

echo ${out_f}
