#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

tab_f='/scratch/groups/mrivas/ukbb24983/phenotypedata/2007183/40831/download/ukb40831.tab'
out_f='/scratch/groups/mrivas/ukbb24983/phenotypedata/extras/assessment_center/ukb40831_field54.tsv.gz'

cols='29-32'

# Field 54: UK Biobank assessment centre
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=54
#
# we use our helper script to identify the tab file as well as the column idx
#
# $ ukbb-query_find_table_by_field_id.sh 54

if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

{
    echo "#FID IID f.54.0.0 f.54.0.1 f.54.0.2 f.54.0.3" | tr ' ' '\t' 
    cut -f1,${cols} ${tab_f} | awk -v OFS='\t' '(NR>1){print $1, $1, $2, $3, $4, $5}' 
} > ${out_f%.gz}

bgzip -l9 -@6 ${out_f%.gz}

