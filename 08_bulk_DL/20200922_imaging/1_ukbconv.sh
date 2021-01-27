#!/bin/bash
set -beEuo pipefail

ml load ukbb-showcase-utils

basket=2005693
table=41413
fids=(20202 20203 20204 20254 20260) # please specify the list of UKB fields here

oak_d=/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/${basket}/${table}/download
cd ${oak_d}

for fid in ${fids[@]} ; do

    bulk_f=/scratch/groups/mrivas/ukbb24983/phenotypedata/${basket}/${table}/bulk/ukb${basket}.${table}.${fid}/ukb${basket}.${table}.${fid}.bulk
    if [ ! -d $(dirname ${bulk_f}) ] ; then mkdir -p $(dirname ${bulk_f}) ; fi
    bulk_oak_f=$(echo ${bulk_f} | sed -e 's%/scratch%/oak/stanford%g')

    if [ ! -f ${bulk_oak_f} ] ; then
        # export the list of individuals with bulk file
        ukbconv ukb${table}.enc_ukb bulk -s${fid}
        mv ukb${table}.bulk ${bulk_f}
        cp ${bulk_f} ${bulk_oak_f}
    else
        # copy the file from OAK to scratch
        cp ${bulk_oak_f} ${bulk_f}
    fi

    # copy the key file
    cp k24983r${table}.key $(dirname ${bulk_f})/.ukbkey
done

