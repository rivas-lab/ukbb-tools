#!/bin/bash
set -beEuo pipefail

plink_gz=$1
out_f=$(dirname $(dirname ${plink_gz}))/wc_l/$(basename $(dirname ${plink_gz}))/$(basename ${plink_gz%.gz}).wc_l.txt

if [ $# -gt 1 ] ; then out_f=$2 ; fi
if [ ! -d $(dirname ${out_f}) ] ; then mkdir -p $(dirname ${out_f}) ; fi

cat_or_zcat () {
    local file=$1
    if [ "${file%.gz}.gz" == "${file}" ] || [ "${file%.bgz}.bgz" == "${file}" ] ; then 
        zcat ${file} 
    elif [ "${file%.zst}.zst" == "${file}" ] ; then 
        zstdcat ${file}
    else
        cat ${file}
    fi
}

{
    echo "#file wc_l"
    echo ${plink_gz} $(cat_or_zcat ${plink_gz} | wc -l) 
} | tr ' ' '\t' > ${out_f}

exit 0
##############
{
    find \
        /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/__to_delete/metal_v2/ \
        /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/metal \
        -type f -name "*.metal.tsv.gz"
} | sort -V > 7e_wc_l_metal.$(date +%Y%m%d-%H%M%S).lst 

6667
7e_wc_l_metal.20201106-091127.lst

sbatch -p mrivas,normal --time=1:0:00 --mem=8000 --nodes=1 --cores=1 \
    --job-name=wc_l --output=logs/wc_l.sh.%A_%a.out --error=logs/wc_l.sh.%A_%a.err \
    --array=1-667 ${parallel_sbatch_no_err_check_sh} 7e_wc_l_metal.sh 7e_wc_l_metal.20201106-091127.lst 10

echo "#file wc_l" | tr ' ' '\t' > /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/metal.wc_l.tsv
find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l/metal -type f \
    | sort -V | parallel -k "cat {} | egrep -v '^#'" >> /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/metal.wc_l.tsv

echo "#file wc_l" | tr ' ' '\t' > /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/__to_delete/metal_v2.wc_l.tsv
find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/wc_l/__to_delete/metal_v2 -type f \
    | sort -V | parallel -k "cat {} | egrep -v '^#'" >> /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE/__to_delete/metal_v2.wc_l.tsv


