#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.0"
NUM_POS_ARGS="1"

source "${SRCDIR}/0_functions.sh"

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
	
	Usage: $PROGNAME [options] output_dir
	  output_dir      The output directory
	
	Options:
	  --cores      (-t)  Number of CPU cores
	
	Default configurations:
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
cores=2
GBE_ID=HC382
n_batch=100
genotype_name=array-combined
out_prefix=ukb24983_v2_hg19
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
        '--GBE_ID' )
            GBE_ID=$2 ; shift 2 ;
            ;;
        '--nbatch' )
            nbatch=$2 ; shift 2 ;
            ;;
        '--genotype_name' )
            genotype_name=$2 ; shift 2 ;
            ;;
        '--out_prefix' )
            out_prefix=$2 ; shift 2 ;
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

out_dir="${params[0]}"

plink_out="${out_dir}/${out_prefix}.${GBE_ID}"
glm_suffix=$(get_plink_suffix ${GBE_ID})

############################################################

# check the log & results files

seq ${n_batch} | tr ' ' '\n' | while read batch_idx ; do
    log_f=${plink_out}.batch${batch_idx}.glm.log
    res_f=${plink_out}.batch${batch_idx}.${glm_suffix}

    if [ ! -s ${log_f} ] || [ ! -s ${res_f} ] ; then
        echo ${batch_idx}
    fi
done

