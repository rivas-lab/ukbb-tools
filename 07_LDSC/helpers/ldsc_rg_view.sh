#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="2.0.0"
NUM_POS_ARGS="0"

############################################################
# functions
############################################################

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Show the tabulated lines for LDSC rg log file
	
	Usage: $PROGNAME [options]
	
	Options:
	  --file (-f)  LDSC rg log file
	  --list (-l)  List of LDSC rg log files
EOF
}

show_LDSC () {
    local LDSC_rg_log=$1
    cat ${LDSC_rg_log} \
    | grep -A2 'Summary of Genetic Correlation Results' \
    | awk 'NR>1' \
    | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'
}

show_LDSC_body () {
    show_LDSC $1 | awk 'NR>1'
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
file=__NONE__
list=/dev/stdin
list_mode=TRUE
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
        '-f' | '--file' )
            file=$2 ; shift 2 ;
            list_mode="FALSE"
            ;;
        '-l' | '--list' )
            list=$2 ; shift 2 ;
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

if [ ${#params[@]} -gt 0 ]; then
    file="${params[0]}"
    list_mode="FALSE"
fi

############################################################

if [ "${list_mode}" == "TRUE" ] ; then
    tmp_list=${tmp_dir}/list.txt
    cat ${list} > ${tmp_list}
    show_LDSC $(cat ${tmp_list} | awk 'NR==1') | awk 'NR==1'
    cat ${tmp_list} | while read f ; do
        ! show_LDSC ${f} | awk 'NR>1'
    done
else
    show_LDSC ${file}
fi
