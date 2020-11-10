#!/bin/bash
set -beEuo pipefail

in_f=$1
pop=$(basename $(dirname ${in_f}))
pheno_name=$(basename ${in_f} | awk -v FS='.' '{print $2}')

cores=3
mem=12000

if [ $# -gt 1 ] ; then cores=$2 ; fi
if [ $# -gt 2 ] ; then mem=$3   ; fi

# in_f=/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/others/ukb24983_exomeOQFE.HC1295.glm.logistic.hybrid.gz
# pop='others'
# pheno_name='HC1295'

out_f=${in_f%.gz}

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf ${tmp_dir} ; }
trap handler_exit EXIT

############################################################
# functions
############################################################

source /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/extras/20201026_exome_gwas_parallel/0_functions.sh

show_remaining_var_lst () {
    local in_f=$1
    local pfile=$2
    local nf=$3
    
    cat_or_zcat_prefix ${in_f} | awk -v FS='\t' -v nf=${nf} '(NF==nf && NR>1){print $3}' \
    | comm --nocheck-order -23 <(cat_or_zcat_prefix ${pfile}.pvar | awk -v FS='\t' '(NR>1){print $3}') /dev/stdin
}

get_nf () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo 13
    else
        echo 14
    fi
}


############################################################
# parameters
############################################################

plink2_version=20201020
pfile=/scratch/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE
# pfile=/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb24983_exomeOQFE
pheno=__AUTO__
keep=__AUTO__
covar_names=__AUTO__
covar_names_add=''
covar=/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe
AUTO_keep=/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983___POP__.phe

pheno_colname=PHENO1

############################################################
# main
############################################################

if [ "${keep}" == "__AUTO__" ] ;        then keep=$(         echo ${AUTO_keep}         | sed -e "s/__POP__/${pop}/g" ) ; fi
if [ "${covar_names}" == "__AUTO__" ] ; then covar_names=$( get_covar_names ${pop} ),${covar_names_add} ; fi

load_plink2 ${plink2_version}
ml load R/3.6 gcc

plink2_glm_remaining_wrapper () {
    local in_f=$1
    local nf=$2    
    shift 2
    
    show_remaining_var_lst ${in_f} ${pfile} ${nf} \
    | plink2 \
    --memory ${mem} \
    --threads ${cores} \
    --pfile ${pfile} $([ -s "${pfile}.pvar.zst" ] && echo "vzs" || echo "") \
    --chr 1-22,X,Y \
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

############################################################


out_dir=${tmp_dir}
if [ ! -d ${out_dir} ] ; then mkdir -p ${out_dir} ; fi

tmp_in_f=${tmp_dir}/$(basename ${in_f%.gz})

! zcat ${in_f} | awk -v FS='\t' -v nf=$(get_nf ${pheno_name}) '(NF==nf)' | uniq > ${tmp_in_f}

tmp_in_f_wc_l=$(cat ${tmp_in_f} | cut -f3 | uniq | wc -l )

if [ ${tmp_in_f_wc_l} -eq 17777817 ] || [ ${tmp_in_f_wc_l} -eq 17777951 ] ; then
    Rscript 8_patch_fix_uniq.R ${tmp_in_f} ${out_f} \
    $([ -s "${pfile}.pvar.zst" ] && echo "${pfile}.pvar.zst" || echo "${pfile}.pvar")
else

    plink_out="${out_dir}/patch"

    plink2_glm_remaining_wrapper ${tmp_in_f} $(get_nf ${pheno_name}) \
    --out ${plink_out} \
    --pheno $([ "${pheno}" == "__AUTO__" ] && get_phe_file ${pheno_name} || echo "${pheno}")

    Rscript 8_patch_fix_sub.R ${tmp_in_f} ${plink_out}.${pheno_colname}.$(get_plink_suffix ${pheno_name}).zst ${out_f} \
    $([ -s "${pfile}.pvar.zst" ] && echo "${pfile}.pvar.zst" || echo "${pfile}.pvar")
fi

bgzip -f -l9 -@${cores} ${out_f}

echo "[$(date +%Y%m%d-%H%M%S)] Patch-fix finished. wc -l: $(zcat ${out_f}.gz | wc -l)"
