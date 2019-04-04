#!/bin/bash
set -beEuo pipefail

PROGNAME=$(basename $(readlink -f $0))
VERSION="2.0"
NUM_POS_ARGS=1

# parse
usage () {
    echo "$PROGNAME (version $VERSION)"
    echo "  Extract phenotype from the master phenotype file"
    echo "  Usage: $PROGNAME [--covar (-c) ] [--covaronly] [--masterphe (-m) master.phe ] GBE_ID"
}

error_option_req_arg (){
   echo "$PROGNAME: option requires an argument -- $1" >&2 ; exit 1
}

FLAG_covar=0
FLAG_covaronly=0
master_phe="$OAK/dev-ukbb-tools/phewas/resources/master.phe"
declare -a params=()
for OPT in "$@" ; do
    case "$OPT" in 
        '-h' | '--help' )
            usage >&2 ; exit 1 ; 
            ;;
        '-v' | '--version' )
            echo $VERSION ; exit 0 ;
            ;;
        '-c' | '--covar' )
            FLAG_covar=1 ; shift 1
            ;;
        '--covaronly' )
            FLAG_covaronly=1 ; shift 1
            ;;
        '-m' | '--masterphe' )
            if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then error_option_req_arg $1 ; fi
            master_phe=$2 ; shift 2
            ;;
        '--'|'-' )
            shift 1 ; params+=( "$@" ) ; break
            ;;
        -*)
            echo "$PROGNAME: illegal option -- '$(echo $1 | sed 's/^-*//')'" 1>&2 ; exit 1
            ;;
        *)
            if [[ ! -z "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
                params+=( "$1" ) ; shift 1
            fi
            ;;
    esac
done

# special case. --covaronly does not require GBE_ID
if [ ${FLAG_covaronly} -eq 1 ] ; then NUM_POS_ARGS=0 ; fi

if [ ${#params[@]} -lt ${NUM_POS_ARGS} ]; then
    echo "$PROGNAME: ${NUM_POS_ARGS} positional arguments are required" >&2
    usage >&2 ; exit 1 ; 
fi


# find the column index of phenotype of interest
if [ ${FLAG_covaronly} -eq 1 ] ; then phe_col_id="" ; else
    GBE_ID=${params[0]}
    phe_col_id=",$( cat $master_phe | awk 'NR==1' | tr "\t" "\n" | egrep -n "${GBE_ID}$" | awk -v FS=':' '{print $1}' )"
fi

# extract the relevant columns
if [ ${FLAG_covaronly} -eq 1 ] || [ ${FLAG_covar} -eq 1 ] ; then
    cut_first_n_fields="45"
else
    cut_first_n_fields="2"
fi

cat $master_phe | cut -f"1-${cut_first_n_fields}${phe_col_id}"

