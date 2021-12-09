#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="1.1.0"
NUM_POS_ARGS="2"

source $(dirname ${SRCDIR})/functions.sh
source "${SRCDIR}/0_functions.sh"


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
	Run gwas

	Usage: $PROGNAME [options] batch_idx output_dir
	  batch_idx       The array job index
	  output_dir      The output directory

	Options:
	  --cores      (-t)  Number of CPU cores
	  --mem        (-m)  The memory amount (MB)

	Default configurations:
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
cores=2
mem=7000
overwrite=FALSE
QT_ALL=FALSE
GBE_ID=HC382
plink2_version=20210701
pop=white_british
n_batch=100
genotype_name=array-combined
out_prefix=ukb24983_v2_hg19
covar_names=__AUTO__
covar_names_add=''
covar=/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe
master_phe=/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/master.20201002.sorted.phe.zst
pheno_col_nums=1994-3562
pheno_colname=PHENO1
pheno=__AUTO__
keep=__AUTO__
pfile=__AUTO__
one_array=__AUTO__
both_arrays=__AUTO__
AUTO_keep=/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983___POP__.phe
AUTO_pfile=/scratch/groups/mrivas/ukbb24983/__GENOTYPE_NAME__/pgen/ukb24983_cal_hla_cnv
AUTO_one_array=/oak/stanford/groups/mrivas/ukbb24983/__GENOTYPE_NAME__/pgen/one_array_variants.txt
AUTO_both_arrays=/oak/stanford/groups/mrivas/ukbb24983/__GENOTYPE_NAME__/pgen/both_arrays_variants.txt
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
        '-t' | '--cores' | '--nCores' )
            cores=$2 ; shift 2 ;
            ;;
        '-m' | '--mem' | '--memory' )
            mem=$2 ; shift 2 ;
            ;;
        '--overwrite' )
            overwrite="TRUE" ; shift 1 ;
            ;;
        '--QT_ALL' )
            QT_ALL="TRUE" ; shift 1 ;
            ;;
        '--GBE_ID' )
            GBE_ID=$2 ; shift 2 ;
            ;;
        '--n_batch' )
            n_batch=$2 ; shift 2 ;
            ;;
        '--plink2_version' )
            plink2_version=$2 ; shift 2 ;
            ;;
        '--pop' )
            pop=$2 ; shift 2 ;
            ;;
        '--genotype_name' )
            genotype_name=$2 ; shift 2 ;
            ;;
        '--out_prefix' )
            out_prefix=$2 ; shift 2 ;
            ;;
        '--covar_names' )
            covar_names=$2 ; shift 2 ;
            ;;
        '--covar_names_add' )
            covar_names_add=$2 ; shift 2 ;
            ;;
        '--covar' )
            covar=$2 ; shift 2 ;
            ;;
        '--master_phe' )
            master_phe=$2 ; shift 2 ;
            ;;
        '--pheno' )
            pheno=$2 ; shift 2 ;
            ;;
        '--pheno_col_nums' )
            pheno_col_nums=$2 ; shift 2 ;
            ;;
        '--keep' )
            keep=$2 ; shift 2 ;
            ;;
        '--pfile' )
            pfile=$2 ; shift 2 ;
            ;;
        '--one_array' )
            one_array=$2 ; shift 2 ;
            ;;
        '--both_arrays' )
            both_arrays=$2 ; shift 2 ;
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

batch_idx="${params[0]}"
out_dir="${params[1]}"

if [ "${pheno}" == "__AUTO__" ] ;       then pheno=$( get_phe_file ${GBE_ID} ) ; fi
if [ "${keep}" == "__AUTO__" ] ;        then keep=$(         echo ${AUTO_keep}         | sed -e "s/__POP__/${pop}/g" ) ; fi
if [ "${pfile}" == "__AUTO__" ] ;       then pfile=$(        echo ${AUTO_pfile}        | sed -e "s/__GENOTYPE_NAME__/${genotype_name}/g" ) ; fi
if [ "${one_array}" == "__AUTO__" ] ;   then one_array=$(    echo ${AUTO_one_array}    | sed -e "s/__GENOTYPE_NAME__/${genotype_name}/g" ) ; fi
if [ "${both_arrays}" == "__AUTO__" ] ; then both_arrays=$(  echo ${AUTO_both_arrays}  | sed -e "s/__GENOTYPE_NAME__/${genotype_name}/g" ) ; fi

n_batch_one_array=$(compute_n_batch_one_array ${n_batch} ${pfile} ${one_array})

if [ "${covar_names}" == "__AUTO__" ] ; then covar_names=$( get_covar_names ${pop} ${batch_idx} ${n_batch_one_array} ),${covar_names_add} ; fi

############################################################

load_plink2 ${plink2_version}
if [ ! -d ${out_dir} ] ; then mkdir -p ${out_dir} ; fi
if [ ! -d ${out_dir}/logs ] ; then mkdir -p ${out_dir}/logs ; fi

plink2_glm_wrapper () {

    show_var_list ${batch_idx} ${n_batch} ${n_batch_one_array} ${one_array} ${both_arrays} \
    | plink2 \
    --memory ${mem} \
    --threads ${cores} \
    --pfile ${pfile} $([ -s "${pfile}.pvar.zst" ] && echo "vzs" || echo "") \
    --chr 1-22,X,XY,Y,MT \
    --covar ${covar} \
    --covar-name $( echo ${covar_names} | tr ',' ' ' ) \
    --extract /dev/stdin \
    --glm zs skip-invalid-pheno firth-fallback cc-residualize hide-covar omit-ref no-x-sex \
    --keep ${keep} \
    --covar-variance-standardize \
    --pheno-quantile-normalize \
    --vif 100000000 \
    $@
}

if [ "${QT_ALL}" == "FALSE" ] ; then

    plink_out="${out_dir}/${out_prefix}.batch${batch_idx}.${GBE_ID}"
    glm_suffix=$(get_plink_suffix ${GBE_ID})
    log_f=$(dirname ${plink_out})/logs/$(basename ${plink_out}).glm.log

    if [ "${overwrite}" == "TRUE" ] || [ ! -s ${log_f} -a ! -s ${plink_out}.${glm_suffix} ] ; then
        plink2_glm_wrapper \
        --out ${plink_out} \
        --pheno ${pheno}

        for ext in ${glm_suffix} ${glm_suffix}.zst ${glm_suffix}.gz ; do
            if [ -f ${plink_out}.${pheno_colname}.${ext} ] ; then
                mv ${plink_out}.${pheno_colname}.${ext} ${plink_out}.${ext}
            fi
        done
        mv ${plink_out}.log ${log_f}

    fi
else
    # QT_ALL mode
    #  apply GWAS scan for multiple quantitative phenotypes
    #  we specify the list of phenotypes in master phenotype file using plink2's --pheno-col-nums
    plink_out="${out_dir}/${out_prefix}.batch${batch_idx}"
    log_f=$(dirname ${plink_out})/logs/$(basename ${plink_out}).QT_ALL.glm.log

    if [ "${overwrite}" == "TRUE" ] || [ ! -s ${log_f} ] ; then
        plink2_glm_wrapper \
        --out ${plink_out} \
        --pheno $(show_scracth_if_exists ${master_phe}) \
        --pheno-col-nums ${pheno_col_nums}

        mv ${plink_out}.log ${log_f}
    fi
fi

