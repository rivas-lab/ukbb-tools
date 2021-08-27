#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="2.0.2"
NUM_POS_ARGS="2"

# source "${SRCDIR}/ldsc_misc.sh"

############################################################
# update log
############################################################
# version 2.0.2 (2020/6/26)
#   With a minor update in the custom pre-processing script, we support the output file from
#   our Metal wrapper script
#
# version 2.0.1 (2020/6/25)
#   We now truncate the p-value at 1e-300
#   https://github.com/bulik/ldsc/issues/144

ml load ldsc
ldscore=${TWINSUK_oak}
merge_alleles=${LDSC_hm3_list}
############################################################
# functions
############################################################

show_default_helper () {
    cat ${SRCNAME} | grep -n Default | tail -n+3 | awk -v FS=':' '{print $1}' | tr "\n" "\t"
}

show_default () {
    cat ${SRCNAME} \
        | tail -n+$(show_default_helper | awk -v FS='\t' '{print $1+1}') \
        | head  -n$(show_default_helper | awk -v FS='\t' '{print $2-$1-1}')
}

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Run ldsc

	Usage: $PROGNAME [options] input_file output_file
	  input_file      The input file
	  output_file     The output file [.sumstats.gz,.log]

	Options:
	  --scratch       Use ldscore in /scratch space
	  --ldscore       Specify the LD score file
	  --merge_alleles Specify the merge allele file (see --merge)
	  --merge         Merge with Hap-map v3 SNP list (you can change the SNP list with --merge_alleles option)

	Note:
	  Please don't use --merge for the summary statistics for the array data (you will lose a lot of variants).

	Default configurations (please use the options above to modify them):
	  ldscore=${ldscore}
	  merge_alleles=${merge_alleles}
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

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
# parser start
############################################################
## == Default parameters (start) == ##
merge="FALSE"
## == Default parameters (end) == ##

declare -a params=()
for OPT in "$@" ; do
    case "$OPT" in
        '-h' | '--help' )
            usage >&2 ; exit 0 ;
            ;;
        '-v' | '--version' )
            echo $VERSION ; exit 0 ;
            ;;
        '--scratch' )
            ldscore=${TWINSUK_scratch} ; shift 1 ;
            ;;
        '--ldscore' )
            ldscore=$2 ; shift 2 ;
            ;;
        '--merge_alleles' )
            merge_alleles=$2 ; shift 2 ;
            ;;
        '--merge' )
            merge="TRUE" ; shift 1 ;
            ;;
        '--'|'-' )
            shift 1 ; params+=( "$@" ) ; break
            ;;
        -*)
            echo "$PROGNAME: illegal option -- '$(echo $1 | sed 's/^-*//')'" 1>&2 ; exit 1
            ;;
        *)
            if [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-+ ]]; then
                params+=( "$1" ) ; shift 1
            fi
            ;;
    esac
done

if [ ${#params[@]} -lt ${NUM_POS_ARGS} ]; then
    echo "${PROGNAME}: ${NUM_POS_ARGS} positional arguments are required" >&2
    usage >&2 ; exit 1 ;
fi

input_file=$(readlink -f "${params[0]}")
output_file=$(readlink -f "${params[1]}")

############################################################

tmp_intermediate_file=${tmp_dir}/$(basename $input_file).input.tsv

echo "Applying a custom pre-processing R script ..."
# echo "${SRCDIR}/make_ldsc_input_file_v2.R /dev/stdout ${input_file} ${ldscore}"

Rscript ${SRCDIR}/make_ldsc_input_file_v2.R /dev/stdout ${input_file} ${ldscore} \
| sed -e 's/[0-9].[0-9][0-9]*[eE]-[1-9][0-9][0-9][0-9][0-9]/1.0e-300/' \
| sed -e 's/[0-9].[0-9][0-9]*[eE]-[1-9][0-9][0-9][0-9]/1.0e-300/' \
| sed -e 's/[0-9].[0-9][0-9]*[eE]-[3-9][0-9][0-9]/1.0e-300/' > ${tmp_intermediate_file}
# truncate the small p-value at 1e-300
# https://github.com/bulik/ldsc/issues/144

echo "Running LDSC munge_sumstats.py ..."

if [ "${merge}" == "TRUE" ] ; then
    munge_sumstats.py \
    --sumstats ${tmp_intermediate_file} \
    --N-col OBS_CT --a1 A1 --a2 A2 --snp ID --signed-sumstats BETA,0 \
    --merge-alleles ${merge_alleles} \
    --out ${output_file%.sumstats.gz}
else
    munge_sumstats.py \
    --sumstats ${tmp_intermediate_file} \
    --N-col OBS_CT --a1 A1 --a2 A2 --snp ID --signed-sumstats BETA,0 \
    --out ${output_file%.sumstats.gz}
fi

echo "output: ${output_file%.sumstats.gz}.sumstats.gz"
