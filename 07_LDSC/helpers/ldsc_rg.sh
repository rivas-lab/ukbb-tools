#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="2.0.0"
NUM_POS_ARGS="3"

# source "${SRCDIR}/ldsc_misc.sh"

ml load ldsc
ldscore=${TWINSUK_oak}
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
	
	Usage: $PROGNAME [options] input_f_1 input_f_2 output_file
	  input_f_1       The 1st input file
	  input_f_2       The 2nd input file
	  output_file     The output file [.log]
	
	Options:
	  --scratch        Use ldscore in /scratch space
	  --ldscore
	
	Default configurations (please use the options above to modify them):
	  ldscore=${ldscore}
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

in1=$(readlink -f "${params[0]}")
in2=$(readlink -f "${params[1]}")
out=$(readlink -f "${params[2]}")

############################################################

tmp_out=${tmp_dir}/$(basename ${out%.log})
out_d=$(dirname ${out})

if [ ! -d ${out_d} ] ; then mkdir -p ${out_d} ; fi

ldsc.py --rg ${in1},${in2} --ref-ld-chr ${ldscore} --w-ld-chr ${ldscore} --out ${tmp_out}

# copy the results to the final output location.
# replace the temp output file path in the log file
cat ${tmp_out}.log | sed -e "s%${tmp_dir}%${out_d}%g" > ${out%.log}.log

echo "Results are written to:  ${out%.log}.log"

