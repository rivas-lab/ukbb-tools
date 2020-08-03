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

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Convert a summary statistics into a PLINK format.

	Usage: $PROGNAME [options] input_file
	  input_file      The input file

	Options:
	  --logit   (-l)  Specify it as a logistic regression.

	Default configurations:
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
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

show_header () {
    local file=$1
    ! cat_or_zcat $file | head -n1
    # cat_or_zcat $file | egrep '^#'
}

get_col_idx () {
    local file=$1
    local key=$2
    show_header $file | sed -e "s/^#//g" | tr "\t" "\n" | awk -v key=$key '($0 == key){print NR}'
}

show_header_linear () {
    echo "#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P ERRCODE" | tr ' ' '\t'
}

show_header_logistic_hybrid () {
    echo "#CHROM POS ID REF ALT A1 FIRTH? TEST OBS_CT OR LOG(OR)_SE Z_STAT P ERRCODE" | tr ' ' '\t'
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
logit=FALSE
col_key_CHROM=CHROM
col_key_POS=POS
col_key_ID=ID
col_key_REF=REF
col_key_ALT=ALT
col_key_A1=A1
col_key_OBS_CT=OBS_CT
col_key_BETA=BETA
col_key_OR=OR
col_key_SE=SE
col_key_P=P
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
        '-l' | '--logit' )
            logit="TRUE" ; shift 1 ;
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
# CHROM POS ID REF ALT A1        TEST OBS_CT BETA SE       T_STAT P ERRCODE
# CHROM POS ID REF ALT A1 FIRTH? TEST OBS_CT OR LOG(OR)_SE Z_STAT P ERRCODE
# CHROM POS ID REF ALT A1             OBS_CT BETA SE              P Direction       HetISq  HetChiSq        HetDf   HetPVal

col_CHROM=$(   get_col_idx $in_file "$col_key_CHROM")
col_POS=$(     get_col_idx $in_file "$col_key_POS")
col_ID=$(      get_col_idx $in_file "$col_key_ID")
col_REF=$(     get_col_idx $in_file "$col_key_REF")
col_ALT=$(     get_col_idx $in_file "$col_key_ALT")
col_A1=$(      get_col_idx $in_file "$col_key_A1")

col_OBS_CT=$(  get_col_idx $in_file "$col_key_OBS_CT")
col_BETA=$(    get_col_idx $in_file "$col_key_BETA")
col_OR=$(      get_col_idx $in_file "$col_key_OR")
col_SE=$(      get_col_idx $in_file "$col_key_SE")
col_P=$(       get_col_idx $in_file "$col_key_P")

if [ "${logit}" == "FALSE" ] ; then
    # linear regression

    show_header_linear

    cat_or_zcat ${in_file} \
    | egrep -v '^#' \
    | awk -v OFS='\t' -v FS='\t' \
        -v col_CHROM=${col_CHROM}   -v col_POS=${col_POS}   -v col_ID=${col_ID} \
        -v col_REF=${col_REF}       -v col_ALT=${col_ALT}   -v col_A1=${col_A1} \
        -v col_OBS_CT=${col_OBS_CT} -v col_BETA=${col_BETA} \
        -v col_SE=${col_SE}         -v col_P=${col_P} \
        -v val_TEST="ADD"           -v val_ERRCODE='.' \
        '{print $col_CHROM, $col_POS, $col_ID, $col_REF, $col_ALT, $col_A1, val_TEST, $col_OBS_CT, $col_BETA, $col_SE, 0, $col_P, val_ERRCODE}'

elif [ "${col_OR}" == "" ] ; then
    # logistic regression, but we compute OR as exp(BETA)

    show_header_logistic_hybrid

    cat_or_zcat ${in_file} \
    | egrep -v '^#' \
    | awk -v OFS='\t' -v FS='\t' \
        -v col_CHROM=${col_CHROM}   -v col_POS=${col_POS}   -v col_ID=${col_ID} \
        -v col_REF=${col_REF}       -v col_ALT=${col_ALT}   -v col_A1=${col_A1} \
        -v col_OBS_CT=${col_OBS_CT} -v col_BETA=${col_BETA} \
        -v col_SE=${col_SE}         -v col_P=${col_P} \
        -v val_TEST="ADD"           -v val_Firth="N"        -v val_ERRCODE="." \
        '{print $col_CHROM, $col_POS, $col_ID, $col_REF, $col_ALT, $col_A1, val_Firth, val_TEST, $col_OBS_CT, exp($col_BETA), $col_SE, 0, $col_P, val_ERRCODE}'
else
    # logistic regression

    show_header_logistic_hybrid

    cat_or_zcat ${in_file} \
    | egrep -v '^#' \
    | awk -v OFS='\t' -v FS='\t' \
        -v col_CHROM=${col_CHROM}   -v col_POS=${col_POS}   -v col_ID=${col_ID} \
        -v col_REF=${col_REF}       -v col_ALT=${col_ALT}   -v col_A1=${col_A1} \
        -v col_OBS_CT=${col_OBS_CT} -v col_OR=${col_OR} \
        -v col_SE=${col_SE}         -v col_P=${col_P} \
        -v val_TEST="ADD"           -v val_Firth="N"        -v val_ERRCODE="." \
        '{print $col_CHROM, $col_POS, $col_ID, $col_REF, $col_ALT, $col_A1, val_Firth, val_TEST, $col_OBS_CT, $col_OR, $col_SE, 0, $col_P, val_ERRCODE}'
fi
