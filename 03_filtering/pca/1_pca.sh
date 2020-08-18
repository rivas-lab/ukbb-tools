#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.2.0"
NUM_POS_ARGS="1"

# source "${SRCDIR}/pca_misc.sh"

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

# var_list_to_bed () {
#     local pvar=$1
#     local one_array=$2
#
#     ml load R/3.6 gcc
#
#     Rscript /dev/stdin ${pvar} ${one_array} << EOF
#     suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
#     args <- commandArgs(trailingOnly=TRUE)
#     pvar <-args[1]
#     one_array <- args[2]    
#     fread(cmd=paste('zstdcat', pvar), colClasses=c('#CHROM'='character')) %>%
#     rename ('CHROM' = '#CHROM') %>%
#     filter(ID %in% (fread(one_array, head=F) %>% pull())) %>%
#     mutate(POSe = POS + length(REF)) %>%
#     select(CHROM, POS, POSe, ID) %>%
#     fwrite('/dev/stdout', sep='\t', col.names=F)    
# EOF
# }

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Run pca
	
	Usage: $PROGNAME [options] out_prefix
	  out_prefix   The prefix of the output files
	
	Options:
	  --keep       (-k)  The set of individuals used in the analysis
	  --pfile            The genotype file in PLINK 2.0 pgen/pvar.zst/psam format
	  --nCores     (-t)  Number of CPU cores
	  --memory     (-m)  The memory amount
	
	Default configurations
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

############################################################
# tmp dir
############################################################
# tmp_dir_root="$LOCAL_SCRATCH"
# if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
# tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# # echo "tmp_dir = $tmp_dir" >&2
# handler_exit () { rm -rf $tmp_dir ; }
# trap handler_exit EXIT

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
nCores=4
memory=30000
keep=""
pfile=/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19
one_array=/oak/stanford/groups/mrivas/ukbb24983/sqc/one_array_variants.txt
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
        '-m' | '--memory' )
            memory=$2 ; shift 2 ;
            ;;
        '--keep' )
            keep=$2 ; shift 2 ;
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

out_prefix="${params[0]}"
# shift 1

############################################################

ml load zstd plink2/20200727 

if [ ! -f $(dirname ${out_prefix}) ] ; then mkdir -p $(dirname ${out_prefix}) ; fi

plink_mem=$( perl -e "use List::Util qw[min max]; print( max( int(${memory} * 0.8), ${memory} - 10000 ))" )
plink_common_opts=" --memory ${plink_mem} --threads ${nCores}"

if [ ! -s ${out_prefix}.prune.log ] && 
   [ ! -s ${out_prefix}.prune.in ] && 
   [ ! -s ${out_prefix}.prune.out ] ; then
#  We apply LD pruning for the variant that passed the following criteria:
#    - autosomal (--chr) biallelic (--max-alleles) QC-passed (--var-filter)
#      common (MAF 5%, --maf, --max-maf) SNPs (--snps-only)
#    - missing rate is at most 10% (--geno)
#    - not significant in HWE test (--hwe)
#    - not in the MHC region (6:25477797-36448354, --exclude)
#      https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37.p13
#      https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37.p13
#    - not in the one array variant set (--exclude)
#  If the keep file is not empty, we focus on the individuals specified in the file (--keep)
        
    plink2 ${plink_common_opts} \
    --pfile ${pfile} $([ -f "${pfile}.pvar.zst" ] && echo "vzs" || echo "") \
    --chr 1-22 --max-alleles 2 --var-filter --snps-only just-acgt \
    --maf 0.05 --max-maf 0.95 \
    --geno 0.1 --hwe 1e-10 midp \
    --rm-dup exclude-all \
    --exclude <( zstdcat ${pfile}.pvar.zst \
      | awk '($1 == 6 && 25477797 <= $2 && $2 < 36448354){print $3}' \
      | cat /dev/stdin ${one_array} \
      ) \
    $([ "${keep}" != "" ] && echo "--keep ${keep}" || echo "" ) \
    --indep-pairwise 50 5 .5 \
    --out ${out_prefix}

    mv ${out_prefix}.log ${out_prefix}.prune.log
fi

# Apply PCA
if [ ! -s ${out_prefix}.eigenvec.log ] && 
   [ ! -s ${out_prefix}.eigenval ] && 
   [ ! -s ${out_prefix}.eigenvec ] && 
   [ ! -s ${out_prefix}.eigenvec.allele.zst ] ; then

    plink2 ${plink_common_opts} \
    --pfile ${pfile} $([ -f "${pfile}.pvar.zst" ] && echo "vzs" || echo "") \
    --extract ${out_prefix}.prune.in \
    $([ "${keep}" == "" ] && echo "" || echo "--keep ${keep}") \
    --pca 40 allele-wts approx vzs vcols=chrom,pos,ref,alt1,alt,ax \
    --seed 24983 \
    --out ${out_prefix}

    mv ${out_prefix}.log ${out_prefix}.eigenvec.log
fi

exit 0
