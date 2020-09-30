#!/bin/bash
set -beEuo pipefail

ml load ukbb-showcase-utils

basket=2005693
table=41413
fids=(20204 20254) # 20202 20203 20260

cd /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/${basket}/${table}/download

#for fid in 20202 20203 20260 ; do
#for fid in 20204 20254 ; do
for fid in ${fids[@]} ; do

    bulk_f=/scratch/groups/mrivas/ukbb24983/phenotypedata/${basket}/${table}/bulk/ukb${basket}.${table}.${fid}/ukb${basket}.${table}.${fid}.bulk
    if [ ! -d $(dirname ${bulk_f}) ] ; then mkdir -p $(dirname ${bulk_f}) ; fi

    # export the list of individuals with bulk file
    ukbconv ukb${table}.enc_ukb bulk -s${fid}
    mv ukb${table}.bulk ${bulk_f}

    # copy the key file
    cp k24983r${table}.key $(dirname ${bulk_f})/.ukbkey
done

