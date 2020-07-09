#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="1"

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

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] || [ "${file%.bgz}.bgz" == "${file}" ] ; then
        zcat ${file}
    elif [ "${file%.zst}.zst" == "${file}" ] ; then
        zstdcat ${file}
    else
        cat ${file}
    fi
}

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Find variant in UKB dataset

	Usage: $PROGNAME [options] CHROM POS

	Options:
      --ID              variant ID (not supported yet)
	  --assembly (-a)   assembly (hg19 or hg38)

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
ID=__AUTO__
assembly=hg19
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
        '--ID' )
            ID=$2 ; shift 2 ;
            ;;
        '-a' | '--assembly' )
            assembly=$2 ; shift 2 ;
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

############################################################

hg38_files=(
/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/ukb24983_exome.pvar.zst
)

hg19_files=(
/oak/stanford/groups/mrivas/ukbb24983/array-combined/pgen/ukb24983_cal_hla_cnv.pvar.zst
/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen/ukb24983_cal_hla_cnv_imp.pvar.zst
)

hg19_imp="/oak/stanford/groups/mrivas/ukbb24983/imp/pgen/ukb24983_imp_chrCHROM_v3.pvar.zst"

if [ ${assembly} == "hg38" ] ; then
    for pvar_f in ${hg38_files[@]} ; do
        echo $pvar_f
        cat_or_zcat $pvar_f | awk -v CHROM=$chrom -v POS=$pos 'NR==1 || ($1 == CHROM && $2 == POS)'
    done
fi

if [ ${assembly} == "hg19" ] ; then
    for pvar_f in ${hg19_files[@]} ; do
        echo $pvar_f
        cat_or_zcat $pvar_f | awk -v CHROM=$chrom -v POS=$pos 'NR==1 || ($1 == CHROM && $2 == POS)'
    done
fi
