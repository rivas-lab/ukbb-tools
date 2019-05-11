#!/bin/bash
set -beEuo pipefail

PROGNAME=$(basename $(readlink -f $0))
VERSION="0.3"

_default_population="white_british"
_default_prefix="ukb24983_v2_hg19"
_default_variant_type="genotyped"

usage () {
    echo "$PROGNAME (version $VERSION)"
    echo "  Check the GWAS sumstats (exists? number of lines and number of non-NA lines)"
    echo "  currently support GWAS on array or exome only"
    echo "  Usage: $PROGNAME [population (default: ${_default_population})] [prefix (default: ${_default_prefix})] [variant_type (default: ${_default_variant_type}]"
}

# read common func
source $(dirname $(readlink -f $0))/04_gwas_misc.sh

# cmdarg parse
if [ $# -lt 1 ] ; then population=${_default_population} ; else population=$1 ;  fi
if [ $# -lt 2 ] ; then prefix=${_default_prefix} ; else prefix=$2 ; fi
if [ $# -lt 3 ] ; then variant_type=${_default_variant_type} ; else variant_type=$3 ; fi

# create a temp directory
tmp_dir="$(mktemp -p /dev/shm/ -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)" 
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

# dump list of phe_files to tmp
tmp_phe_files=$tmp_dir/phe_files.lst
cat ../05_gbe/phenotype_info.tsv | awk -v min_N=10 'NR > 1 && $7 >= min_N' | egrep -v MED  | awk '{print $NF}' > $tmp_phe_files

tmp_parallel_wrapper=$tmp_dir/parallel_wrapper.sh
echo "#!/bin/bash" > ${tmp_parallel_wrapper}
echo "source $(dirname $(readlink -f $0))/04_gwas_misc.sh" >> ${tmp_parallel_wrapper}
echo 'check_sumstats_wrapper $@' >> ${tmp_parallel_wrapper}

echo "#idx phe_file sumstats_file n_lines n_non-NA_lines" | tr " " "\t"
parallel -j 6 --keep "bash $tmp_parallel_wrapper $tmp_phe_files {} $population $prefix $variant_type" ::: $(seq 1 $(cat  $tmp_phe_files | wc -l))

