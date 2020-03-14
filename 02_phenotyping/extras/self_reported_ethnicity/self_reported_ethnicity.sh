#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.0"
NUM_POS_ARGS="1"

source "${SRCDIR}/self_reported_ethnicity_misc.sh"

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
	Extract self-reported ethnicity (field 21000)
	
	Usage: $PROGNAME [options] basket_id table_id
	  basket_id         The basket ID
      table_id          The table ID
	
	Options:
	  --out_dir         The output file directory
	
	Default configurations:
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
out_dir=/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/self_reported_ethnicity/misc
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
        '--out_dir' )
            out_dir=$2 ; shift 2 ;
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

basket_id="${params[0]}"
table_id="${params[1]}"

############################################################

out_tsv="${out_dir}/ukb${basket_id}_ukb${table_id}_f21000.tsv"
out_phe="$(dirname $(dirname ${out_tsv}))/phe/$(basename ${out_tsv%.tsv}.phe)"

tab_file=$(find_tab_file ${basket_id} ${table_id})
tab_col_file="${tab_file}.columns"

tab_file_copy=$(find_scratch_copy ${tab_file})

cols=$(find_col_idx_by_field_id ${tab_file} 21000)

############################################################

if [ ! -f ${out_tsv} ] ; then
    extract_cols ${tab_file_copy} ${cols} > ${out_tsv}
fi
if [ ! -f ${out_phe} ] ; then
    Rscript ${SRCDIR}/self_reported_ethnicity_phe.R ${out_tsv} ${out_phe}
fi
