#!/bin/bash
set -beEuo pipefail

ml load ukbb-showcase-utils

cd /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/41413/download

for fid in 20202 20203 20260 ; do

    bulk_f=/scratch/groups/mrivas/ukbb24983/phenotypedata/2005693/41413/bulk/ukb2005693.41413.${fid}/ukb2005693.41413.${fid}.bulk

    # export the list of individuals with bulk file
    ukbconv ukb41413.enc_ukb bulk -s${fid}
    mv ukb41413.bulk ${bulk_f}

    # copy the key file
    cp k24983r41413.key $(dirname ${bulk_f})/.ukbkey
done

