#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.0"
NUM_POS_ARGS="3"

# read common func
source "$(dirname ${SRCNAME})/09_liftOver_misc.sh"

############################################################
# functions
############################################################
show_default_helper () {
    cat ${SRCNAME} | grep -n 'Default_parameters_(' | tail -n+2 \
    | awk -v FS=':' '{print $1}' | tr "\n" "\t" 
}

show_default () {
    cat ${SRCNAME} \
        | tail -n+$(show_default_helper | awk -v FS='\t' '{print $1+1}') \
        | head  -n$(show_default_helper | awk -v FS='\t' '{print $2-$1-1}')
}

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Apply UCSC liftOver
	  Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
	
	Usage: $PROGNAME [options] in_file out_mapped out_unmapped
	  in_file       The plink2 pgen/pvar.zst/psam file.
	  out_mapped    The phenotype file
	  out_unmapped  The name of the phenotype. We assume the phenotype is stored with the same column name
	
	Options:
	  --threads  (-t) Number of CPU cores
	  --src_genome    The genome build for the input file
	  --dst_genome    The genome build for the output file
	
	Default configurations (please use the options above to modify them):
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
## == Default_parameters_(start) == ##
threads=4
src_genome="hg19"
dst_genome="hg38"
## == Default_parameters_(end) == ##

declare -a params=()
for OPT in "$@" ; do
    case "$OPT" in 
        '-h' | '--help' )
            usage >&2 ; exit 0 ; 
            ;;
        '-v' | '--version' )
            echo $VERSION ; exit 0 ;
            ;;
        '-t' | '--threads' )
            threads=$2 ; shift 2 ;
            ;;
        '--src_genome' )
            src_genome=$2 ; shift 2 ;
            ;;
        '--dst_genome' )
            dst_genome=$2 ; shift 2 ;
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

in_file="${params[0]}"
out_mapped="${params[1]}"
out_unmapped="${params[2]}"

############################################################
liftOverWrapper ${in_file} ${src_genome} ${dst_genome} ${out_mapped} ${out_unmapped} ${threads}
