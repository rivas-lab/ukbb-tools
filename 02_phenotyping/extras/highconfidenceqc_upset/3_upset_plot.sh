#!/bin/bash
set -beEuo pipefail

# ml load R/3.6 gcc inkscape

GBE_ID=$1

if [ $# -gt 1 ] ; then
    pop=$2
else
    pop="none"
fi
# GBE_ID='HC294'
# pop='white_british'

data_d="/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc_upset"
HC_idx_long_f="${data_d}/ukb37855_ukb40831_icd.annot.tsv.gz"

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

if [ "${pop}" == "" ] ||
   [ "${pop}" == "NONE" ] ||
   [ "${pop}" == "None" ] ||
   [ "${pop}" == "none" ] ; then

    keep_f='none'
    out_pdf="${data_d}/upset_plot/all/${GBE_ID}.pdf"
else

    keep_f="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_${pop}.phe"
    out_pdf="${data_d}/upset_plot/${pop}/${GBE_ID}.${pop}.pdf"
fi

if [ ! -d $(dirname ${out_pdf}) ] ; then mkdir -p $(dirname ${out_pdf}) ; fi


if [ -s "${out_pdf}" ] && [ -s "${out_pdf%.pdf}.png" ] ; then
    exit 0
fi

# extract the relevant
in_tbl=${tmp_dir}/in_tbl.tsv
tabix -h ${HC_idx_long_f} ${GBE_ID} > ${in_tbl}

# calling Rscript in Yosuke's singularity image.
# you should be able to run it with your own R environment by simply droppiing run-simg.sh
run-simg.sh Rscript /dev/stdin ${in_tbl} ${keep_f} ${out_pdf} << EOF
suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
library(UpSetR)

args <- commandArgs(trailingOnly=TRUE)

in_tbl  <- args[1]
keep_f  <- args[2]
out_pdf <- args[3]

####################
# function
####################

upsetplot_wrapper <- function(upset_labels, HC_idx_long, text_scaling_factor=1){
    upset(
        fromList(as.list(setNames(
            upset_labels %>% lapply(function(l){
            # extract the list of individuals
                HC_idx_long %>%
                filter(upset_label == l) %>%
                pull(IID) %>% unique()
            }),
            upset_labels
        ))),
        order.by = "freq", show.numbers = "yes",
        nsets = 20, nintersects = 40,
        text.scale = text_scaling_factor * c(1.5, 1.2, 1.5, 1.2, 1, .6),
#         number.angles = 300,
#         point.size = 2, line.size = .5,
#         mb.ratio = c(0.6, 0.4),
        mainbar.y.label = "Number of case individuals",
        sets.x.label = "# cases per data source"
    )
}

####################
# main
####################

HC_idx_long <- fread(in_tbl) %>%
rename('GBE_ID'='#GBE_ID') %>%
select(-time, -array) %>%
mutate(
    upset_label = if_else(coding == 6, 'Self-reported', paste('ICD-10', val))
)

if(keep_f != 'none'){
# read a keep file and focus on the subset of specified individuals

    keep_f %>%
    fread(sep='\t', col.names=c('FID', 'IID'), colClasses = 'character') %>%
    pull(FID) -> keep_l

    HC_idx_long %>%
    filter(IID %in% keep_l) -> HC_idx_long

}

# enumerate the list of phenotype sources
HC_idx_long %>%
count(upset_label) %>%
arrange(-n) %>%
pull(upset_label) -> upset_labels

cairo_pdf(out_pdf, height = 6, width = 8, family = "Helvetica")
upsetplot_wrapper(upset_labels, HC_idx_long)
dev.off()

EOF

# convert pdf to png
inkscape --export-filename=${out_pdf%.pdf}.png ${out_pdf}

echo ${out_pdf}
echo ${out_pdf%.pdf}.png
