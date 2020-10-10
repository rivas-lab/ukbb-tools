#!/bin/bash
set -beEuo pipefail


idx_f=$1


{
    cat /scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003/cnv.20201006.unfinished.pvar | egrep '^#'

    cat $idx_f | awk 'length($0)>0'  | while read idx ; do
        cat /scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003/input_20201006/split.cnv.20201006.unfinished.pvar.body.${idx}.pvar | egrep -v '#'
    done
} > /scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003/cnv.20201007.unfinished.pvar

