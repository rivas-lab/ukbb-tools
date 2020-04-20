#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.1"
NUM_POS_ARGS="2"

source "${SRCDIR}/18_metal_misc.sh"

flipcheck_sh="$(dirname ${SRCDIR})/09_liftOver/flipcheck.sh"

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
	Run run_metal for the specified set of summary statistics and apply flipfix.
	
	Usage: $PROGNAME [options] in_file_list outfile_prefix
	  in_file_list      A file that has a list of input files for METAL
	  metal_out_prefix  The prefix of output files from METAL
	
	Options:
	  --flipcheck_sh     The location of flip check script
	  --nCores     (-t)  Number of CPU cores
	  --assembly    The genome build for the input file (option for flipcheck)
	  --ref_fa      The reference genome sequence. (option for flipcheck)
	
	Note:
	  We assume the input file has the following columns:
	    - OR, CHROM, POS, ID, A1, REF, BETA, P, SE, OBS_CT
	  The output files will be
	    - <outfile_prefix>.metal.tsv.gz : METAL output file
	    - <outfile_prefix>.metal.info.txt : log file
      This script internally calls flipcheck.sh. Please check 09_liftOver for more info.
	
	Default configurations:
	  flipcheck_sh=${flipcheck_sh}
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
# trap handler_exit EXIT

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
nCores=4
ref_fa="AUTO"
assembly="hg19"
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
        '-t' | '--nCores' )
            nCores=$2 ; shift 2 ;
            ;;
        '--flipcheck_sh' )
            flipcheck_sh=$2 ; shift 2 ;
            ;;
        '--ref_fa' )
            ref_fa=$2 ; shift 2 ;
            ;;
        '--assembly' )
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

input_file="${params[0]}"
out_file="${params[1]}"

############################################################

ml load metal

tmp_out=${tmp_dir}/$(basename ${out_file})
master_file="${tmp_dir}/metal.masterfile"

show_master_file ${input_file} ${tmp_out} > ${master_file}

cd ${tmp_dir}
metal ${master_file}
cd -

extract_loci_for_files ${input_file} ${nCores} | bgzip -l9 -@ ${nCores} > ${tmp_out}.loci.gz

Rscript ${SRCDIR}/metal_post_processing_step1.R \
${tmp_out}.loci.gz ${tmp_out}1.tbl ${tmp_out}.metal.tsv

bash ${flipcheck_sh} --ref_fa ${ref_fa} --assembly ${assembly} ${tmp_out}.metal.tsv \
| bgzip -l 9 -@ ${nCores} > ${tmp_out}.metal.check.tsv.gz

Rscript ${SRCDIR}/metal_post_processing_step2.R \
${tmp_out}.metal.check.tsv.gz ${tmp_out}.metal.fixed.tsv

bgzip -l 9 -@ ${nCores} ${tmp_out}.metal.fixed.tsv

if [ ! -d $(dirname ${out_file}) ] ; then mkdir -p $(dirname ${out_file}) ; fi

cp ${tmp_out}.metal.fixed.tsv.gz ${out_file}.metal.tsv.gz
cat ${tmp_out}1.tbl.info ${master_file} > ${out_file}.metal.info.txt
echo "the results are written in: ${out_file%.gz}.metal.tsv.gz"

