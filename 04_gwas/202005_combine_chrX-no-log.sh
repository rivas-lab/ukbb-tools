#!/bin/bash
set -beEuo pipefail

# A helper script to combine the autosomal and chrX GWAS runs.
#
# Yosuke Tanigawa (2020/6/2)
#
# In the 2020/05 GWAS refresh, we pushed two sets of jobs, corresponding to the
# autosomal runs and chrX (more precisely, chrX, XY, Y, and MT).
#
# This script performs the followings:
#  1. move the autosome run into new file
#  2. concatenate autosome runs and chrX runs file (sumstats and log files) and save it under the original filename for automal run
#  3. bgzip the combined file
#  4. delete the sym links corresponding to the chrX run.
#
# As a result, the symlinks used to point to the autosomal results point to the one for the combined one
#
# Usage:
#  bash 202005_combine_chrX.sh <autosomal_sumstats_symlink> <chrX_sumstats_symlink>
# 
#    Please pass the symlinks
#
# Example:
#  bash 202005_combine_chrX.sh /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/white_british/ukb24983_v2_hg19.INI50.array-combined.glm.linear.gz /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current/white_british/ukb24983_v2_hg19.INI50_X.array-combined.glm.linear.gz
#
#
# Dependencies:
#  bgzip and R
#

#ml load snpnet_yt # this is Yosuke's R env
#ml load htslib

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)

dataset="array-combined"

# read args
#if [ $# -lt 2 ] ; then 
#    echo "usage: ${PROGNAME} chrAUTO.sumstats chrX.sumstats" >&1 
#    exit 1
#fi
#chrAUTO_sym=$1
#chrX_sym=$2

if [ $# -lt 2 ] ; then 
    echo "usage: ${PROGNAME} pop GBE_ID" >&1 
    exit 1
fi
pop=$1
GBE_ID=$2

GBE_CAT=$(echo $GBE_ID | sed -e "s/[0-9]//g")
if [ ${GBE_CAT} == "INI" ] || [ ${GBE_CAT} == "QT_FC" ] ; then
    glm_suffix="linear"
else
    glm_suffix="logistic.hybrid"
fi

data_dir="/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current"
chrAUTO_sym="${data_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}.${dataset}.glm.${glm_suffix}.gz"
chrX_sym="${data_dir}/${pop}/ukb24983_v2_hg19.${GBE_ID}_X.${dataset}.glm.${glm_suffix}.gz"

# check dir
if [ "$(dirname $chrAUTO_sym)" != "$(dirname $chrX_sym)" ] ; then
    echo "two sumstats (symlinks) are not in the same dir!" >&1; exit 1
fi
symlink_d=$(dirname $chrAUTO_sym)
chrAUTO=$(readlink -f $chrAUTO_sym)
chrX=$(readlink -f $chrX_sym)
if [ "$(dirname $chrAUTO)" != "$(dirname $chrX)" ] ; then
    echo "two sumstats (real path) are not in the same dir!" >&1; exit 1
fi

# check if symlinks are specified in cmdargs
if [ "$chrAUTO" == "$chrAUTO_sym" ] ; then
    echo "Please pass the symlink (not the full path) to this script. $chrAUTO_sym" >&1; exit 1
fi
#if [ "$chrX" == "$chrX_sym" ] ; then
#    echo "Please pass the symlink (not the full path) to this script. $chrX_sym" >&1; exit 1
#fi

# check prefix
# example: ukb24983_v2_hg19.INI50
prefix=$(basename ${chrAUTO_sym%.gz} | sed -e "s/.logistic.hybrid//g" | sed -e "s/.linear//g" | sed -e "s/.glm//g" | sed -e "s/.${dataset}//g" )
prefixX=$(  basename ${chrX_sym%.gz} | sed -e "s/.logistic.hybrid//g" | sed -e "s/.linear//g" | sed -e "s/.glm//g" | sed -e "s/.${dataset}//g")

if [ "${prefix}_X" != "${prefixX}" ] ; then
    echo "The prefix of the file names are not identical: ${prefix}_X ${prefixX}" >&1; exit 1
fi

chrAUTO_log_sym=$symlink_d/logs/${prefix}.${dataset}.log
#chrX_log_sym=$symlink_d/logs/${prefix}_X.${dataset}.log

chrAUTO_log=$(readlink -f $chrAUTO_log_sym)
#chrX_log=$(readlink -f $chrX_log_sym)

#if [ "$(dirname $chrAUTO_log)" != "$(dirname $chrX_log)" ] ; then
#    echo "two log files (real path) are not in the same dir!" >&1; exit 1
#fi

if [ "$(basename $chrAUTO)" != "$(basename $chrAUTO_sym)" ] ; then
    echo "basename mismatch! $chrAUTO" >&1; exit 1
fi
if [ "$(basename $chrX)" != "$(basename $chrX_sym)" ] ; then
    echo "basename mismatch! $chrX" >&1; exit 1
fi
#if [ "$(basename $chrAUTO_log)" != "$(basename $chrAUTO_log_sym)" ] ; then
#    echo "basename mismatch! $chrAUTO_log" >&1; exit 1
#fi
#if [ "$(basename $chrX_log)" != "$(basename $chrX_log_sym)" ] ; then
#    echo "basename mismatch! $chrX_log" >&1; exit 1
#fi

#chrAUTO_log_new=$(dirname ${chrAUTO_log})/${prefix}_autosomes.${dataset}.log
chrAUTO_new=$(dirname ${chrAUTO})/${prefix}_autosomes.$(basename ${chrAUTO} | sed -e "s/${prefix}.//g")

#mv $chrAUTO_log $chrAUTO_log_new
mv $chrAUTO     $chrAUTO_new

# combine the log files
#cat $chrAUTO_log_new $chrX_log > $chrAUTO_log

# combine the sumstats files
Rscript /dev/stdin $chrAUTO_new $chrX ${chrAUTO%.gz} << EOF
args <- commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))
in_f_auto <- args[1]
in_f_X    <- args[2]
out_f <- args[3]

df <- bind_rows(
    fread(in_f_auto, colClasses=c('#CHROM'='character', 'POS'='integer', 'ID'='character', 'P'='character')),
    fread(in_f_X,    colClasses=c('#CHROM'='character', 'POS'='integer', 'ID'='character', 'P'='character'))
) %>%
rename('CHROM'='#CHROM') %>%
mutate(CHROM =  factor(CHROM, levels = c(1:22, 'X', 'XY', 'Y', 'MT'))) %>%
arrange(CHROM, POS) %>%
mutate(CHROM =  as.character(CHROM)) %>%
rename('#CHROM' = 'CHROM') %>%
fwrite(out_f, sep='\t', na = "NA", quote=F)
EOF

# bgzip
bgzip -l9 -f ${chrAUTO%.gz}

# remove the sym links to chrX file
#rm $chrX_log_sym
rm $chrX_sym

