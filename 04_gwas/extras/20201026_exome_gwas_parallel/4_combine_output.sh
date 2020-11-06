#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.0"
NUM_POS_ARGS="1"

source "${SRCDIR}/0_functions.sh"

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
	Run gwas
	
	Usage: $PROGNAME template_f combined_f
	  output_dir      The output directory
	
	Options:
	  --cores      (-t)  Number of CPU cores
	  --n_batch
	
	Default configurations:
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
cores=1
n_batch=100
log=FALSE
check=TRUE
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
        '-t' | '--cores' | '--nCores' )
            cores=$2 ; shift 2 ;
            ;;
        '--n_batch' )
            n_batch=$2 ; shift 2 ;
            ;;
        '--log' )
            log="TRUE" ; shift 1 ;
            ;;
        '--skip_check' )
            skip_check="FALSE" ; shift 1 ;
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

template_f="${params[0]}"
combined_f="${params[1]}"

############################################################

if [ "${check}" == "TRUE" ] ; then
    combine_check_files ${template_f} ${n_batch} > ${tmp_dir}/check_files.txt
    if [ $(cat ${tmp_dir}/check_files.txt | wc -l  ) -gt 0 ] ; then
        echo "We found missing files" 2>&1
        cat ${tmp_dir}/check_files.txt
        exit 1
    fi
fi

if [ "${log}" == "TRUE" ] ; then
    combine_log_files   ${template_f} ${n_batch} ${combined_f} ${cores}
else
    combine_plink_files ${template_f} ${n_batch} ${combined_f} ${cores}
fi

echo ${combined_f}

