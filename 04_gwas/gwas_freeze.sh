#!/bin/bash
set -beEuo pipefail

__default_pop="white_british"
__default_p_val_thr='1e-3'
__default_out_dir_root="/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/freeze/"

show_usage () {
    echo "$0 check_file prefix variant_type freeze_v [pop (default: ${__default_pop})] [p_val_threshold (default: ${__default_p_val_thr})] [out_dir_root (default: ${__default_out_dir_root})]"
    echo ""
    echo "example: $0 check_array_gwas.XXXXXXX.out ukb24983_v2_hg19 genotyped 20190512 white_british"
}

if [ $# -lt 4 ] ; then show_usage >&2 ; exit 1 ; fi

check_file=$(readlink -f $1)
prefix=$2
variant_type=$3
freeze_v=$4
if [ $# -gt 4 ] ; then pop=$5 ;          else pop=${__default_pop} ; fi
if [ $# -gt 5 ] ; then p_val_thr=$6 ;    else p_val_thr=${__default_p_val_thr} ; fi
if [ $# -gt 6 ] ; then out_dir_root=$7 ; else out_dir_root=${__default_out_dir_root} ; fi

#####################################

show_gwas_src="$(dirname $(readlink -f $0))/show_gwas.sh"

p_val_thr_str=$(echo ${p_val_thr} | sed -e 's/-//g')
master_prefix="${prefix}.${pop}.${variant_type}.glm"
master_file="${master_prefix}.${freeze_v}.${p_val_thr_str}.tsv"

out_dir=${out_dir_root}/${freeze_v}/${pop}

if [ ! -d ${out_dir} ] ; then mkdir -p ${out_dir} ; fi

cd ${out_dir}

cp -a ${check_file} .
echo "tar"
cat ${check_file} | grep -v NA | cut -f3 | tar --transform 's/.*\///g' --group root --owner root -cvf "${master_prefix}.${freeze_v}.tar" -T /dev/stdin
echo "extract as unsorted"
bash ${show_gwas_src} ${check_file} ${p_val_thr} > ${master_file%.tsv}.unsorted.tsv
echo "header"
cat ${master_file%.tsv}.unsorted.tsv | awk 'NR==1' > ${master_file}
echo "sort"
cat ${master_file%.tsv}.unsorted.tsv | awk 'NR>1' | sort -S 60G --parallel=4 -k1V,1 -k2n,2 -k4V,4 >> ${master_file}
echo "bgzip"
bgzip --compress-level 9 -f ${master_file}
echo "tabix"
tabix ${master_file}.gz -s 1 -b 2 -e 2
echo "symlinks"
ln -s ${master_file}.gz     ${master_prefix}.${p_val_thr_str}.tsv.gz 
ln -s ${master_file}.gz.tbi ${master_prefix}.${p_val_thr_str}.tsv.gz.tbi
rm ${master_file%.tsv}.unsorted.tsv

