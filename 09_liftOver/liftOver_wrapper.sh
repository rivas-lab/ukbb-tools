#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.1"
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
	Apply UCSC liftOver. If the input file contains OR or BETA column, the script will also automatically apply the flipfix.
	  Author: Yosuke Tanigawa (ytanigaw@stanford.edu)
	
	Usage: $PROGNAME [options] in_file out_mapped out_unmapped
	  in_file       The input file. Please see the notes blow about our assumptions on the input file format.
	  out_mapped    The output file of the mapped elements.
	  out_unmapped  The output file of the unmapped elements.
	
	Options:
	  --threads  (-t) Number of CPU cores
	  --src_genome    The genome build for the input file
	  --dst_genome    The genome build for the output file
	
	Notes:
	  - We assume the input file has a header line (first line) that starts with `#`.
	  - The input file needs to have the following columns: CHROM, POS, REF, and ID. 
	  - If effect size column (BETA or OR) is provided, we will automatically apply flipfix. 
	  - You may include additional columns.
	
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
out_mapped="${params[1]}"
out_unmapped="${params[2]}"

############################################################
ml load ucsc-liftover
liftOverWrapper ${in_file} ${src_genome} ${dst_genome} ${out_mapped} ${out_unmapped} ${threads} ${to_bed_field_sep} ${bed_chr_prefix}
