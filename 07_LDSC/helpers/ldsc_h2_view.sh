#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="0"

############################################################
# functions
############################################################

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Show the tabulated lines for LDSC h2 log file

	Usage: $PROGNAME [options]

	Options:
	  --file (-f)  LDSC h2 log file
	  --list (-l)  List of LDSC h2 log files
EOF
}

show_LDSC_h2 () {
    local LDSC_h2_log=$1
    echo "#p h2_obs h2_obs_se lambda_GC mean_chi2 intercept intercept_se ratio ratio_se" | tr ' ' '\t'

    local trait=$(cat ${LDSC_h2_log} | grep 'Reading summary statistics from' | sed -e "s/Reading summary statistics from//g" | awk '{print $1}' )

    cat ${LDSC_h2_log} \
    | grep -A4 'Total Observed scale h2' \
    | sed -e 's/Ratio < 0/Ratio: <0/g' \
    | awk -v FS=':' '{print $2}' \
    | sed -e 's/[()]//g' \
    | sed -e 's/mean chi^2 < 1/mean_chi^2<1/g' \
    | sed -e 's/usually indicates GC correction./usually_indicates_GC_correction/g' \
    | tr '\n' ' ' \
    | awk -v p=${trait} -v OFS='\t' '{print p, $1, $2, $3, $4, $5, $6, $7, $8}'
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
    show_LDSC_h2 $(cat ${tmp_list} | awk 'NR==1') | awk 'NR==1'
    cat ${tmp_list} | while read f ; do
        ! show_LDSC_h2 ${f} | awk 'NR>1'
    done
else
    show_LDSC_h2 ${file}
fi
