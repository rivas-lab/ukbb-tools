#!/bin/bash
set -beEuo pipefail

bulk_file=$1
out_dir=$(readlink -f $2)

# create a temp directory
tmpDir=$(mktemp -d -p $LOCAL_SCRATCH tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX) ; echo "temp_dir = $tmpDir"
handler_exit () { rm -rf $tmpDir ; }
trap handler_exit EXIT

bulk_file_missing=${tmpDir}/$(basename $bulk_file)-$(date +%Y%m%d-%H%M%S)

find ${out_dir} -type f \
| awk -v FS='/' '{print $NF}' \
| awk -v FS='.' '{print $1}' \
| sort -k1 \
| comm -1 -3 /dev/stdin <( cat ${bulk_file} | awk -v OFS="_" '{print $1, $2}'| sort -k1 ) \
| awk -v FS='_' -v OFS='_' '{print $1 " " $2, $3, $4}' > ${bulk_file_missing}

n_lines=$(cat ${bulk_file_missing} | wc -l)
n_batches=$( perl -e "print(int((${n_lines} + 999)/1000))" )

echo "Downloading ${n_lines} files (${bulk_file_missing})"

parallel -j 4 $(dirname $(readlink -f $0))/$(basename $0 .sh)_sub.sh ${bulk_file_missing} {} ${out_dir} :::  $(seq 1 ${n_batches}) 

