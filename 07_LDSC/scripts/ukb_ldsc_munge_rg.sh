#!/bin/bash
set -beEuo pipefail

show_usage_and_exit () {
        echo "usage: $0 in_PLINK_1 in_PLINK_2 dataset_name [out_dir_root]" >&1 ; exit 1
}

_default_out_dir_root="$OAK/dev-ukbb-tools/ldsc/"
repo_dir=$(dirname $(dirname $(readlink -f $0)))

infer_in_type () {
    sumstats=$1
    sumstats_base=$(basename $sumstats .gz)
    if [ $(   echo ${sumstats_base} | rev | cut -c -16 | rev ) == ".logistic.hybrid" ] ; then
	echo "logistic"
    elif [ $( echo ${sumstats_base} | rev | cut -c -7  | rev ) == ".linear" ] ; then
	echo "linear"
    fi
}

get_GBE_name () {
    sumstats=$1
    sumstats_base=$(basename $sumstats .gz)
    echo ${sumstats_base} | awk -v FS='.' '{print $2}'
}

set_rg_file_base () {
    sumstats1=$1
    sumstats2=$2
    sumstats_base1=$(basename $sumstats1 .gz)
    sumstats_base2=$(basename $sumstats2 .gz)
    prefix=$(echo ${sumstats_base1} | awk -v FS='.' -v OFS='.' '{print $1, $3}')
    echo "${prefix}.$(get_GBE_name $sumstats1)_$(get_GBE_name $sumstats2).log" 
}

get_munged_file () {
    sumstats=$1
    dataset_name=$2
    out_dir_root=$3
    echo "${out_dir_root}/munged/${dataset_name}/$(basename $sumstats | awk -v FS='.' -v OFS='.' '{print $1, $2, $3}').sumstats.gz"
}

get_rg_file () {
    sumstats1=$1
    sumstats2=$2
    dataset_name=$3
    out_dir_root=$4
    echo "${out_dir_root}/rg/${dataset_name}/$(set_rg_file_base ${sumstats1} ${sumstats2})"
}

run_munged () {
    in_file=$1
    dataset_name=$2
    out_dir_root=$3
    out_file=$(get_munged_file $in_file $dataset_name $out_dir_root)
    if [ ! -f ${out_file}.log ] ; then	
       if [ ! -d $(dirname ${out_file}) ] ; then mkdir -p $(dirname ${out_file}) ; fi
        echo bash ${repo_dir}/helpers/ukb_ldsc_munge_wrapper.sh \
	    ${out_file} $(get_GBE_name $in_file) $(infer_in_type $in_file) $in_file
        bash ${repo_dir}/helpers/ukb_ldsc_munge_wrapper.sh \
	    ${out_file} $(get_GBE_name $in_file) $(infer_in_type $in_file) $in_file
    fi
    echo ${out_file}.log
}

run_rg () {
    in_file1=$1
    in_file2=$2
    dataset_name=$3
    out_dir_root=$4
    munged1=$(get_munged_file $in_file1 $dataset_name $out_dir_root)
    munged2=$(get_munged_file $in_file2 $dataset_name $out_dir_root)
    out_file=$(get_rg_file $in_file1 $in_file2 $dataset_name $out_dir_root)
    if [ ! -f ${out_file} ] ; then
        if [ ! -d $(dirname ${out_file}) ] ; then mkdir -p $(dirname ${out_file}) ; fi    
        bash ${repo_dir}/helpers/ukb_ldsc_rg_helper.sh \
   	    ${out_file} ${munged1} ${munged2}
    fi
    echo ${out_file}
}


if [ $# -lt 3 ] ; then show_usage_and_exit ; fi
in_sumstats1=$( readlink -f $1 )
in_sumstats2=$( readlink -f $2 )
dataset_name=$3
if [ $# -gt 3 ] ; then out_dir_root=$4 ; else out_dir_root=${_default_out_dir_root} ; fi


run_munged ${in_sumstats1} ${dataset_name} ${out_dir_root} 
run_munged ${in_sumstats2} ${dataset_name} ${out_dir_root} 
run_rg ${in_sumstats1} ${in_sumstats2} ${dataset_name} ${out_dir_root}

