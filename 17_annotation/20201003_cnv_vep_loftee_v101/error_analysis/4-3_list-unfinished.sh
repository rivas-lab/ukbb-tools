#!/bin/bash
set -beEuo pipefail

out=$(basename $0 .sh).$(date +%Y%m%d-%H%M%S).lst

find /scratch/groups/mrivas/ukbb24983/cnv/annotation_20201003/output_vep_20201006 -name "*.tsv" \
    | rev | awk -v FS='.' '{print $3}' | rev | sort  | comm -3 <(seq -w 9851 | sort)  /dev/stdin \
    | sort -nr > ${out}

echo ${out}
