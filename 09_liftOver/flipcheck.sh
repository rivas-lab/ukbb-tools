#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="1.0.0"
NUM_POS_ARGS="1"

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
	Partially apply flipfix and check the flip.
	  Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
	
	Usage: $PROGNAME [options] in_file
	  in_file       The input file. It needs to have the following columns: CHROM, POS, ID, REF, ALT, A1, and one of the following: OR or BETA.
	
	Options:
	  --assembly    The genome build for the input file
	  --ref_fa      The reference genome sequence. It it's specified as AUTO, it will be automatically grab one from ${REF_FA_DIR} (the default behavior).
	
	Default configurations (please use the options above to modify them):
	  to_bed_field_sep=${TO_BED_FIELD_SEP}
	  bed_chr_prefix=${BED_CHR_PREFIX}
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

############################################################
# parser start
############################################################
to_bed_field_sep=${TO_BED_FIELD_SEP}
bed_chr_prefix=${BED_CHR_PREFIX}
## == Default_parameters_(start) == ##
ref_fa="AUTO"
assembly="hg19"
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
        '--ref_fa' )
            ref_fa=$2 ; shift 2 ;
            ;;
        '--assembly' )
            assembly=$2 ; shift 2 ;
            ;;
        '--to_bed_field_sep' )
            to_bed_field_sep=$2 ; shift 2 ;
            ;;
        '--bed_chr_prefix' )
            bed_chr_prefix=$2 ; shift 2 ;
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

############################################################
ml load bedtools

if [ "${ref_fa}" == "AUTO" ] ; then ref_fa=$(get_ref_fa "${assembly}"); fi
check_flip ${in_file} ${ref_fa} ${to_bed_field_sep} ${bed_chr_prefix}
