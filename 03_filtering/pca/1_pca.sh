#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.3.0"
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
	Run PCA with plink2 or smartpca in Eigensoft (--EIG option)
	
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

approx_str () {
    keep=$1
    n_samples=$(cat $keep | egrep -v '^#' | wc -l)
    if [ "${n_samples}" -gt 5000 ] ; then
        echo "approx"
    else
        echo ""
    fi
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
nCores=4
memory=30000
keep=""
EIG="FALSE"
nPCs=40
pfile=/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2_hg19
snp_qc=/oak/stanford/groups/mrivas/ukbb24983/snp/snp_download/ukb_snp_qc.txt
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
        '--nPCs' )
            nPCs=$2 ; shift 2 ;
            ;;
        '--EIG' )
            EIG="TRUE" ; shift 1 ;
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

ml load zstd plink2/20200727 R/3.6 gcc

if [ ${EIG} == "TRUE" ] ; then
    ml load EIG
fi

if [ ! -f $(dirname ${out_prefix}) ] ; then mkdir -p $(dirname ${out_prefix}) ; fi

plink_mem=$( perl -e "use List::Util qw[min max]; print( max( int(${memory} * 0.8), ${memory} - 10000 ))" )
plink_common_opts=" --memory ${plink_mem} --threads ${nCores}"

# Apply PCA

if [ ${EIG} == "TRUE" ] ; then

    if [ ! -s ${out_prefix}.eigensoft.log ] && 
       [ ! -s ${out_prefix}.eigenval ] && 
       [ ! -s ${out_prefix}.eigenvec ] ; then

        bfile=${tmp_dir}/$(basename ${pfile})

        # generate PLINK1.9 BED file

        plink2 ${plink_common_opts} \
        --pfile ${pfile} $([ -f "${pfile}.pvar.zst" ] && echo "vzs" || echo "") \
        --extract <(cat ${snp_qc} | awk '($118 == 1){print $1}' ) \
        $([ "${keep}" == "" ] && echo "" || echo "--keep ${keep}") \
        --make-bed \
        --out ${bfile}

        cp ${bfile}.log ${out_prefix}.bed.log    

        # prepare snp and ind file
        # ref: https://github.com/DReichLab/EIG/tree/master/CONVERTF

        cat ${bfile}.bim | awk '{print $2, $1, "0.0", $4, $6, $5}' > ${bfile}.snp
        cat ${bfile}.fam | awk -v FID_IID_sep=':' -v label=$([ "${keep}" == "" ] && echo "sample" || basename "${keep}") '{print $1 FID_IID_sep $2, "U", label }' > ${bfile}.ind

        # run PCA
        # ref: https://github.com/DReichLab/EIG/tree/master/EIGENSTRAT

        ! smartpca.perl \
        -i ${bfile}.bed \
        -a ${bfile}.snp \
        -b ${bfile}.ind \
        -m 5 \
        -k ${nPCs} \
        -o ${out_prefix}.pca \
        -e ${out_prefix}.evals \
        -l ${out_prefix}.eigensoft.log \
        -p ${out_prefix}.plot

        # post-processing
        { 
            echo "#FID IID $(seq ${nPCs} | awk -v prefix=PC '{print prefix $1}' | tr '\n' ' ' | rev | cut -c2- | rev ) label" |  tr ' ' '\t'

            cat ${out_prefix}.pca.evec | sed -e 's/^\s\+//g' | egrep -v '^#' | sed -e 's/\s\+/\t/g' | tr ':' '\t'
        } > ${out_prefix}.eigenvec

        cat ${out_prefix}.pca | awk '1 < NR' | awk -v nPCs=${nPCs} 'NR <= nPCs' > ${out_prefix}.eigenval

    fi

else

    if [ ! -s ${out_prefix}.eigenvec.log ] && 
       [ ! -s ${out_prefix}.eigenval ] && 
       [ ! -s ${out_prefix}.eigenvec ] && 
       [ ! -s ${out_prefix}.eigenvec.allele.zst ] ; then

        plink2 ${plink_common_opts} \
        --pfile ${pfile} $([ -f "${pfile}.pvar.zst" ] && echo "vzs" || echo "") \
        --extract <(cat ${snp_qc} | awk '($118 == 1){print $1}' ) \
        $([ "${keep}" == "" ] && echo "" || echo "--keep ${keep}") \
        --pca ${nPCs} allele-wts $([ "${keep}" == "" ] && echo "" || $(approx_str ${keep}) ) vzs vcols=chrom,pos,ref,alt1,alt,ax \
        --seed 24983 \
        --out ${out_prefix}

        mv ${out_prefix}.log ${out_prefix}.eigenvec.log
    fi

fi

# plot the first 2 PCs

if [ ! -s ${out_prefix}.eigenvec.PC1.PC2.png ] ; then

    Rscript /dev/stdin ${out_prefix}.eigenvec ${out_prefix}.eigenvec.PC1.PC2.png << EOF
        suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
        args <- commandArgs(trailingOnly=TRUE)

        evec_f <- args[1]
        evec_p <- args[2]

        evec <- fread(evec_f, colClasses=c('#FID'='character', 'IID'='character')) %>%
        rename('FID'='#FID')

        p1 <- evec %>%
        ggplot(aes(x=PC1, y=PC2)) +
        stat_density_2d(aes(fill = ..level..), geom = "polygon") +
        labs(title = file.path(basename(dirname(evec_f)), basename(evec_f))) +
        theme_bw()


        p2 <- evec %>%
        ggplot(aes(x=PC1, y=PC2)) +
        geom_point(alpha=.05) +
        labs(title = file.path(basename(dirname(evec_f)), basename(evec_f))) +
        theme_bw()
        
        ggsave(evec_p, gridExtra::arrangeGrob(p1, p2, ncol=2), width=12, height=6)        
EOF

fi

exit 0
