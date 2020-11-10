#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"


data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"

check_res_out_f=$1

out_f=${PROGNAME%.sh}.$(date +%Y%m%d-%H%M%S).tsv

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
############################################################



{
    echo "#population GBE_ID" | tr ' ' '\t'

    for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' ; do
        if [ -d ${data_d}/${pop}/logs ] ; then
            ! find ${data_d}/${pop}/logs -type f -name "*.log.gz" | grep -v QT_ALL | sort -V | awk -v pop=${pop} -v FS='.' -v OFS='\t' '{print pop, $(NF-3)}'
        fi
    done

} > ${tmp_dir}/finished.tsv

{
    echo "#population GBE_ID" | tr ' ' '\t'

    cat ${check_res_out_f} | awk -v OFS=':' '($3==100){print $1, $2}' | sort \
    | comm -3 /dev/stdin <(cat ${tmp_dir}/finished.tsv | egrep -v '#' | tr '\t' ':' | sort) \
    | tr ':' '\t'

} > ${out_f}

echo ${out_f}
