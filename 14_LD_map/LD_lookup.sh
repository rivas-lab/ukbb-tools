#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="2"

# source "${SRCDIR}/LD_lookup_misc.sh"

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
	Look at LD matrix and show r2
	
	Usage: $PROGNAME [options] CHROM POS
	  CHROM           chromosome
	  POS             position
	
	Options:
	  --r2               The r2 threshold
      --pop              Population
      --ld               The ld file (default: AUTO -- automatically configured)
	
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
r2=0.1
pop=white_british
ld=AUTO
ukb_dir=/oak/stanford/groups/mrivas/ukbb24983
dataset=array_imp_combined_no_cnv
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
        '--r2' )
            r2=$2 ; shift 2 ;
            ;;
        '--pop' )
            pop=$2 ; shift 2 ;
            ;;
        '--ld' )
            ld=$2 ; shift 2 ;
            ;;
        '--ukb_dir' )
            ukb_dir=$2 ; shift 2 ;
            ;;
        '--dataset' )
            dataset=$2 ; shift 2 ;
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

chrom="${params[0]}"
pos="${params[1]}"

############################################################

if [ "${ld}" == "AUTO" ] ; then
    ld="${ukb_dir}/${dataset}/ldmap/ukb24983_cal_hla_imp.${pop}.ld_map.tsv.gz"
fi

tabix -h ${ld} ${chrom}:${pos} \
| awk -v chr=${chrom} -v pos=${pos} -v rsq=${r2} -v OFS='\t' \
'(NR==1){print $0}
 ($1 == chr && $2 == pos && $7 >= rsq){print $1, $2, $3, $4, $5, $6, $7}
 ($4 == chr && $5 == pos && $7 >= rsq){print $4, $5, $6, $1, $2, $3, $7}'
