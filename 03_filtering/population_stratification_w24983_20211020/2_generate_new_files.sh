#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})

ukb_d='/oak/stanford/groups/mrivas/ukbb24983'

if [ ! -s ${ukb_d}/sqc/population_stratification_w24983_20211020/ukb24983_master_sqc.20211020.phe ] ; then
    set -x
#     run-simg.sh is Yosuke's container for his R env, you can replace this with Rscript
    run-simg.sh Rscript \
    add_nonWBsplit_col.R \
    ${ukb_d}/sqc/population_stratification_w24983_20200828/ukb24983_master_sqc.20200828.phe \
    ${ukb_d}/sqc/population_stratification_w24983_20211020/ukb24983_master_sqc.20211020.phe
    set +x
fi

if [ ! -s ${ukb_d}/phenotypedata/master_phe/master.20211020.phe.gz ] ; then
    set -x
    run-simg.sh Rscript \
    add_nonWBsplit_col.R \
    ${ukb_d}/phenotypedata/master_phe/master.20210129.phe.gz \
    ${ukb_d}/phenotypedata/master_phe/master.20211020.phe

    bgzip -l9 -@6 ${ukb_d}/phenotypedata/master_phe/master.20211020.phe
    set +x
fi
