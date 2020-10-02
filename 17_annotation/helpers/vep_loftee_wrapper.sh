#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.1.0"
NUM_POS_ARGS="3"

############################################################
# default values
############################################################
public_d="/oak/stanford/groups/mrivas/public_data"
vep_data="${public_d}/vep/20200912"
loftee_data=$(dirname ${vep_data})/20201002_loftee_data

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
	Run Ensembl's Variant Effect Predictor (VEP) with LOFTEE Plugin

	Usage: $PROGNAME [options] assembly input_vcf_file vep_output_prefix
	  assembly            (GRCh37 or GRCh38)
	  input_vcf_file      The input vcf file
	  vep_output_prefix   The output file prefix

	Options:
	  --memory     (-m)  The memory amount
	  --vep_data         The reference data file for VEP
	  --loftee_data      The reference data directory for loftee plugin
	  --fasta            The reference sequence file

	Default configurations:
	  vep_data=${vep_data}
	  loftee_data=${loftee_data}
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
fasta=__AUTO__
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
        '--vep_data' )
            vep_data=$2 ; shift 2 ;
            ;;
        '--loftee_data' )
            loftee_data=$2 ; shift 2 ;
            ;;
        '--fasta' )
            fasta=$2 ; shift 2 ;
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

assembly="${params[0]}"
vep_in_vcf="${params[1]}"
vep_out="${params[2]}"

############################################################

# load the vep module
# see: https://github.com/rivas-lab/sherlock-modules/tree/master/vep
ml load vep

LoF_d="${loftee_data}/${assembly}"
if   [ "${assembly}" == "GRCh37" ] ; then
    LoF_data_paths="human_ancestor_fa:${LoF_d}/human_ancestor.fa.gz,conservation_file:${LoF_d}/phylocsf_gerp.sql,gerp_file:${LoF_d}/GERP_scores.final.sorted.txt.gz"
elif [ "${assembly}" == "GRCh38" ] ; then
    # this needs to be updated when we test with GRCh38
    LoF_data_paths="human_ancestor_fa:${LoF_d}/human_ancestor.fa.gz,conservation_file:${LoF_d}/phylocsf_gerp.sql,gerp_file:${LoF_d}/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
else
    "error: unsupported assembly (${assembly})" >&2 ; exit 1
fi

if [ "${fasta}" == "__AUTO__" ] ; then
    if   [ "${assembly}" == "GRCh37" ] ; then
        fasta=${vep_data}/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
    elif [ "${assembly}" == "GRCh38" ] ; then
        fasta=${vep_data}/homo_sapiens/101_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
    else
        "error: unsupported assembly (${assembly})" >&2 ; exit 1
    fi
fi

# generate the output directory
if [ ! -d $(dirname ${vep_out}) ] ; then mkdir -p $(dirname ${vep_out}) ; fi

vep \
    --offline --cache \
    --dir_cache ${vep_data} \
    --fasta ${fasta} \
    --allele_number --everything \
    --vcf \
    --assembly ${assembly} -i ${vep_in_vcf} -o ${vep_out} \
    --plugin LoF,loftee_path:/opt/vep/src/loftee-master,${LoF_data_paths}
