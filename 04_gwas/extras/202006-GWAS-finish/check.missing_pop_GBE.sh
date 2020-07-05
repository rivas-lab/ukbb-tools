#!/bin/bash
set -beEuo pipefail

############################################################
# filenames
############################################################

info_file="/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/05_gbe/phenotype_info.tsv"
gwas_current_dir="/oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current"
pops=('white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others')
min_N=100
missing_pop_GBE_f=missing_pop_GBE.minN${min_N}.$(date +%Y%m%d-%H%M%S).tsv

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
# functions
############################################################

get_plink_suffix () {
    local GBE_ID=$1

    GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")

    if [ "${GBE_CAT}" == "QT_FC" ] || [ "${GBE_CAT}" == "INI" ] ; then
        echo glm.linear
    else
        echo glm.logistic.hybrid
    fi
}

show_GBE_IDs () {
    info_file=$1
    N_thr=$2
    cat $info_file | egrep -v '^#' | awk -v FS='\t' -v thr=$N_thr '($7 >= thr && length($1)>0){print $1}'
}

get_sumstats_link () {
    gwas_current_dir=$1
    pop=$2
    GBE_ID=$3

    echo "${gwas_current_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.$(get_plink_suffix ${GBE_ID}).gz"
}

############################################################
# functions
############################################################

# check R env
Rscript /dev/stdin << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
EOF

tmp_missing_symlinks_f=${tmp_dir}/missing.symlinks.tsv

echo "[$(date +%Y%m%d-%H%M%S)] Checking the symlinks in ${gwas_current_dir} ..." >&2

{
    echo "#population GBE_ID"
    for pop in ${pops[@]} ; do
        echo $pop >&2
        show_GBE_IDs ${info_file} ${min_N} | while read GBE_ID ; do
            sumstats_link=$(get_sumstats_link ${gwas_current_dir} ${pop} ${GBE_ID})
            if [ ! -h ${sumstats_link} ] ; then echo $pop $GBE_ID ; fi
        done
    done
} | tr " " "\t" | awk -v FS='\t' 'NF==2' > ${tmp_missing_symlinks_f}

echo "[$(date +%Y%m%d-%H%M%S)] Applying the population-specific N >= ${min_N} filter ..." >&2

Rscript /dev/stdin ${info_file} ${tmp_missing_symlinks_f} ${min_N} ${missing_pop_GBE_f} << EOF
args <- commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
info_f            <- args[1]
check_symlinks_f  <- args[2]
min_N             <- as.integer(args[3])
missing_pop_GBE_f <- args[4]

info_df <- fread(info_f) %>% rename('GBE_ID'='#GBE_ID')

pop_name_map <- data.frame(
    pop=c('N_GBE', 'N_NBW', 'N_AFR', 'N_EAS', 'N_SAS', 'N_SMR', 'N_OTH'),
    population=c('white_british', 'non_british_white', 'african', 's_asian', 'e_asian', 'related', 'others'),
    stringsAsFactors=F
)

check_symlinks_df <- fread(check_symlinks_f) %>% rename('population'='#population')

check_symlinks_df %>%
left_join(
    info_df %>%
    select(GBE_ID, N_GBE, N_NBW, N_AFR, N_EAS, N_SAS, N_SMR, N_OTH) %>%
    gather(pop, N_pop, -GBE_ID) %>%
    left_join(pop_name_map, by='pop') %>%
    select(GBE_ID, population, N_pop),
    by=c('GBE_ID', 'population')
) %>% filter(N_pop >= min_N) %>%
left_join(info_df %>% select(GBE_ID, GBE_NAME), by=c('GBE_ID')) %>%
rename('#population' = 'population') %>%
fwrite(missing_pop_GBE_f, sep='\t', na = "NA", quote=F)
EOF

echo "[$(date +%Y%m%d-%H%M%S)] Results are written in ${missing_pop_GBE_f}" >&2
