#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.2.2"
NUM_POS_ARGS="1"

source "${SRCDIR}/18_metal_misc.sh"

flipcheck_sh="$(dirname ${SRCDIR})/09_liftOver/flipcheck.sh"

############################################################
# update log
############################################################
# version 0.2.1 (2020/6/25)
#   Add OBS_CT column in the Metal output file
#
# version 0.2.1 (2020/6/24)
#   Robust parsing of the Metal output files in post-processing scripts
#
# version 0.2.0 (2020/6/12)
#   Automatic filtering of the NA-lines
#
# version 0.1.2 (2020/4/20)
#   Starting this version, the METAL post-processing script
#   properly handles flipfix for HLA and CNV datasets in UKB

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
	
	Usage: $PROGNAME [options] in_file_1 [in_file_2..n]
	  in_file_1..n     Input files for METAL
	
	Options:
      --out (-o) [REQUIRED] The prefix of output files
	  --flipcheck_sh     The location of flip check script
	  --nCores     (-t)  Number of CPU cores
	  --in_file          List of input files for METAL
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
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
# trap handler_exit EXIT

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
nCores=4
ref_fa="AUTO"
assembly="hg19"
in_file="AUTO"
out="__REQUIRED__"
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
        '-f' | '--in_file' )
            in_file=$2 ; shift 2 ;
            ;;
        '-o' | '--out' )
            out=$2 ; shift 2 ;
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

if [ "${in_file}" == "AUTO" ] && [ ${#params[@]} -lt ${NUM_POS_ARGS} ]; then
    echo "${PROGNAME}: ${NUM_POS_ARGS} positional arguments are required" >&2
    usage >&2 ; exit 1 ; 
fi

if [ "${out}" == "__REQUIRED__" ] ; then
    echo "Output file (--out option) is required" >&2
    usage >&2 ; exit 1
fi

############################################################

ml load metal
ml load R/3.6 gcc/6

# check if we have R env with packages
Rscript /dev/stdin << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
EOF

tmp_infile_dir=${tmp_dir}/inputs
tmp_infile_original=${tmp_dir}/metal.original.input.lst
tmp_infile_list=${tmp_dir}/metal.input.lst
tmp_out=${tmp_dir}/$(basename ${out})
master_file="${tmp_dir}/metal.masterfile"

if [ ! -d ${tmp_infile_dir} ] ; then mkdir -p ${tmp_infile_dir} ; fi

# copy the list of input files
if [ ${in_file} == "AUTO" ] ; then
    for f in ${params[@]} ; do
        if [ "${f}" != "" ] ; then
            echo $f >> ${tmp_infile_original}
        fi
    done
else
    cp ${in_file} ${tmp_infile_original}
fi

cat ${tmp_infile_original} | awk '{print NR, $1}' | while read nr f ; do
    if [ "${f}" != "" ] ; then
        echo "pre-processing $f"
        tmp_f=${tmp_infile_dir}/metal_input.${nr}.gz
        metal_pre_processing ${f} | bgzip -@ ${nCores} > ${tmp_f}
        echo ${tmp_f} >> ${tmp_infile_list}
    fi
done

echo "Generating master file for Metal ..."

# generate METAL master file
show_master_file ${tmp_infile_list} ${tmp_out} > ${master_file}

echo "Applying Metal ..."
cd ${tmp_dir}
metal ${master_file}
cd -

echo "Metal is done. Applying post-processing scripts ..."

echo "Joining the metal output with CHROM, POS, and OBS_CT ..."

# metal_post_processing_step1.R --> join the metal output with CHROM, POS, and OBS_CT columns
Rscript ${SRCDIR}/metal_post_processing_step1.R ${tmp_infile_list} ${tmp_out}1.tbl ${tmp_out}.metal.tsv

echo "Applying flipcheck script to fetch the REF allele from FASTA file ..."

# apply flipcheck script to fetch the REF allele from FASTA file
bash ${flipcheck_sh} --ref_fa ${ref_fa} --assembly ${assembly} ${tmp_out}.metal.tsv | bgzip -@ ${nCores} > ${tmp_out}.metal.check.tsv.gz

echo "Applying flipfix using a custom script ..."
Rscript ${SRCDIR}/metal_post_processing_step2.R ${tmp_out}.metal.check.tsv.gz ${tmp_out}.metal.fixed.tsv

echo "Copying the results ..."
# bgzip
bgzip -l 9 -@ ${nCores} ${tmp_out}.metal.fixed.tsv

# copy the results from tmp_dir to the actual output dir
if [ ! -d $(dirname ${out}) ] ; then mkdir -p $(dirname ${out}) ; fi
cp ${tmp_out}.metal.fixed.tsv.gz ${out%.gz}.metal.tsv.gz
cat <(echo "# == original input files ==") ${tmp_infile_original} <(echo "# == METAL info file ==") ${tmp_out}1.tbl.info <(echo "# == METAL master file ==") ${master_file} <(echo "# == output file ==") <(echo "${out%.gz}.metal.tsv.gz") <(echo "${out%.gz}.metal.info.txt") > ${out%.gz}.metal.info.txt

echo "the results are written in: ${out%.gz}.metal.tsv.gz"
echo "the log file is written to: ${out%.gz}.metal.info.txt"
