#!/bin/bash
set -beEuo pipefail

PROGNAME=$(basename $(readlink -f $0))
VERSION="0.1"

# read common func
source $(dirname $(readlink -f $0))/04_gwas_misc.sh
_default_p_val_thr='1e-4'

usage () {
    echo "$PROGNAME (version $VERSION)"
    echo "  Generate GWAS sumstats file (unsorted)"
    echo "  Usage: $PROGNAME <check_array_gwas.XXXX.out> [p-val threshold (default: ${_default_p_val_thr}]"
}

# cmdarg parse
if [ $# -lt 1 ] ; then usage >&2 ; exit 1 ; fi
in_file=$1
if [ $# -lt 2 ] ; then p_val_thr=${_default_p_val_thr} ; else p_val_thr=$2 ; fi

#in_file='check_array_gwas.42238988.out'
#p_val_thr='1e-4'

# create a temp directory
tmp_dir="$(mktemp -p /dev/shm/ -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)" 
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

tmp_parallel_wrapper=$tmp_dir/parallel_wrapper.sh
echo "#!/bin/bash" > ${tmp_parallel_wrapper}
echo "source $(dirname $(readlink -f $0))/04_gwas_misc.sh" >> ${tmp_parallel_wrapper}
echo 'show_sumstats $@' >> ${tmp_parallel_wrapper}

show_sumstats_header | awk '{print "#" $0}'
cat_or_zcat $in_file | awk 'NR>1' | cut -f 2,3 \
    | while read phe sumstats ; do echo "$(basename ${phe%.phe}) $sumstats" ; done \
    | parallel -j6 --keep "bash $tmp_parallel_wrapper {} ${p_val_thr}" :::: /dev/stdin

