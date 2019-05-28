#!/bin/bash
set -beEuo pipefail

PROGNAME=$(basename $(readlink -f $0))
VERSION="1.0"
NUM_POS_ARGS=1
DEBUG=1 # 0 means true

usage () {
    echo "$PROGNAME (version $VERSION)"
    echo "  Compute column-wise md5check sum for a given table file"
    echo "  Usage: $PROGNAME [--cols (-c) comma-deliminated-col-index ] table.file"
}

error_option_req_arg () {
   echo "$PROGNAME: option requires an argument -- $1" >&2 ; exit 1
}

debug_msg () {
    if [ ${DEBUG} -eq 0 ] ; then echo "[$PROGNAME] [DEBUG] $@" ; fi
}

log_lvl_n () {
    local verbose_level=$1
    local verbose_threshold=$2
    shift 2
    if [ ${verbose_level} -gt ${verbose_threshold} ] || [ ${DEBUG} -eq 0 ] ; then echo "[$PROGNAME] $@" ; fi
}

log_lvl1 () {
    local verbose_level=$1
    shift 1
    log_lvl_n ${verbose_level} 0 $@
}

cat_or_zcat () {
    local in_file=$1
    if [ "${in_file%.gz}.gz" == ${in_file} ] ; then zcat ${in_file} ; else cat ${in_file} ; fi
}

extract_col () {
    local in_file=$1
    local col_idx=$2
    cat_or_zcat ${in_file} | cut -f ${col_idx}
}

extract_col_data () {
    local in_file=$1
    local col_idx=$2
    extract_col ${in_file} ${col_idx} | awk 'NR>1 && length($0) > 0' | grep -v NA
}

get_num_cols () {
    local in_file=$1
    cat_or_zcat ${in_file} | awk 'NR==1' | tr "\t" "\n" | wc -l
}

get_col_name () {
    local in_file=$1
    local col_idx=$2
    extract_col ${in_file} ${col_idx} | awk 'NR==1'
}

compute_col_md5sum () {
    local in_file=$1
    local col_idx=$2
    extract_col ${in_file} ${col_idx} | md5sum | awk '{print $1}'
}

count_non_NAs () {
    local in_file=$1
    local col_idx=$2
    extract_col_data ${in_file} ${col_idx} | wc -l
}

compute_summary () {
    local in_file=$1
    local col_idx=$2
    local n_vals=$3
    if [ $n_vals -gt 10 ] ; then
        extract_col_data ${in_file} ${col_idx} | Rscript \
            -e 'd<-scan("stdin", quiet=TRUE)' \
            -e 'cat(min(d, na.rm=T), median(d, na.rm=T), max(d, na.rm=T), mean(d, na.rm=T), sd(d, na.rm=T), sep="/")'
    else
       extract_col_data ${in_file} ${col_idx} | sort | uniq -c \
           | awk -v OFS=':' -v ORS=',' '{print $2, $1}' | sed -e 's/,$//g'
    fi 
}

show_col_wise_md5sum_per_line () {
    local in_file=$1
    local col_idx=$2
    n_vals=$(extract_col_data ${in_file} ${col_idx} | sort -u | wc -l )
    summary_str=$(compute_summary ${in_file} ${col_idx} ${n_vals}| tr "\n" " " )
    debug_msg "compute_summary = ${n_vals}, summary_str = ${summary_str}" >&2
    if [ "${summary_str}" == "" ] ; then summary_str="-" ; fi
    if [ ${n_vals} -gt 10 ] ; then inferred="continuous" ; else inferred="categorical" ; fi
    echo $(get_col_name ${in_file} ${col_idx}) \
        $(compute_col_md5sum ${in_file} ${col_idx}) \
        $(count_non_NAs ${in_file} ${col_idx}) \
        ${n_vals} \
        ${inferred} \
        ${summary_str}
}

show_col_wise_md5sum () {
    local tmp_in_file=$1
    local n_cols=$2
    local n_arg_cols=$3
    if [ $# -gt 3 ] ; then
        local arg_cols=$4
    fi
    # loop over cols and show the results
    if [ ${n_arg_cols} -eq 0 ] ; then
        for col_idx in $( seq 1 ${n_cols} ) ; do
            echo $(basename ${tmp_in_file}) ${col_idx} $(show_col_wise_md5sum_per_line ${tmp_in_file} ${col_idx})
        done
    else
        for loop_idx in $( seq 1 ${n_arg_cols} ) ; do
            col_idx=$( echo ${arg_cols} | tr "," "\n" | awk -v idx=${loop_idx} 'NR==idx' )
            debug_msg "loop_check ${loop_idx} ${col_idx} ${n_cols}" >&2
            if [ "${col_idx}" -gt 0 ] && [ "${col_idx}" -le "${n_cols}" ]  ; then
                echo $(basename ${tmp_in_file}) ${col_idx} $(show_col_wise_md5sum_per_line ${tmp_in_file} ${loop_idx})
            fi
        done
    fi | tr " " "\t"
}

##########################################################################
# start parse
##########################################################################
verbose_level=0
arg_cols=""
tmp_root=${LOCAL_SCRATCH}
col_start=0
col_end=0
header=1
declare -a params=()
for OPT in "$@" ; do
    case "$OPT" in 
        '-h' | '--help' )
            usage >&2 ; exit 1 ; 
            ;;
        '-v' | '--version' )
            echo $VERSION ; exit 0 ;
            ;;
        '-V' | '--verbose' )
            verbose_level=1 ; shift 1
            ;;
        '-H' | '--header' )
            header=0 ; shift 1
            ;;
        '-c' | '--cols' )
            if [[ $# -lt 2 ]] || [[ "$2" =~ ^-+ ]]; then error_option_req_arg $1 ; fi
            arg_cols=$2 ; shift 2
            ;;
        '-s' | '--col_start' )
            if [[ $# -lt 2 ]] || [[ "$2" =~ ^-+ ]]; then error_option_req_arg $1 ; fi
            col_start=$2 ; shift 2
            ;;
        '-e' | '--col_end' )
            if [[ $# -lt 2 ]] || [[ "$2" =~ ^-+ ]]; then error_option_req_arg $1 ; fi
            col_end=$2 ; shift 2
            ;;
        '-t' | '--tmp_root' )
            if [[ $# -lt 2 ]] || [[ "$2" =~ ^-+ ]]; then error_option_req_arg $1 ; fi
            tmp_root=$2 ; shift 2
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
    exit 1 ; 
fi

in_file=$(readlink -f ${params[0]})
##########################################################################
# end parse
##########################################################################

# create a temp directory
tmp_dir="$(mktemp -d -p ${tmp_root} tmp-${PROGNAME}-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
handler_exit () { rm -rf ${tmp_dir} ; }
trap handler_exit EXIT
log_lvl1 ${verbose_level} "tmp_dir = ${tmp_dir}" >&2

# count number of columns
n_cols=$( cat_or_zcat ${in_file} | awk -v FS='\t' '(NR==1){print NF}'  )
if [ "${arg_cols}" == "" ] ; then 
    n_arg_cols=0
else
    n_arg_cols=$( echo $arg_cols | tr "," "\n" | wc -l )
fi

# col range to cols
debug_msg "range_debug" ${col_start} ${col_end} ${n_arg_cols} >&2
if [ ${col_start} -ne 0 ] && [ ${col_end} -ne 0 ] && [ ${n_arg_cols} -eq 0 ] ; then
    arg_cols=$( seq ${col_start} ${col_end} | tr "\n" "," | sed -e 's/,$//g' )
    (( n_arg_cols = ${col_end} - ${col_start} + 1 ))
fi
log_lvl1 ${verbose_level} "n_cols = ${n_cols} n_arg_cols = ${n_arg_cols}" >&2

# create a tmp copy of the input file
tmp_in_file="${tmp_dir}/$(basename ${in_file})"
if [ ${n_arg_cols} -eq 0 ] ; then 
    cat ${in_file} | awk 'NR==1' > ${tmp_in_file}
    cat ${in_file} | awk 'NR>1' | sort -k1n,1 >> ${tmp_in_file}
else
    cat ${in_file} | awk 'NR==1' | cut -f${arg_cols}  > ${tmp_in_file}
    paste <(cat ${in_file} | cut -f1) <(cat ${in_file} | cut -f ${arg_cols}) | awk 'NR>1' | sort -k1n,1 | cut -f2- >> ${tmp_in_file}
fi
log_lvl1 ${verbose_level} "copy of the input file is made to ${tmp_in_file}" >&2

if [ ${header} -eq 0 ] ; then
    echo "#file col_idx col_name md5sum n_non_NAs n_uniq_vals inferred_data_type summary" | tr " " "\t"
fi

# loop over cols and show the results
show_col_wise_md5sum ${tmp_in_file} ${n_cols} ${n_arg_cols} ${arg_cols}

