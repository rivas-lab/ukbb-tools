#!/bin/bash
set -beEuo pipefail

show_usage_and_exit () {
	echo "usage: $0 in_PLINK_sumstats [out_file_munged] [out_file_h2] [in_type]" >&1 ; exit 1
}

configure_out_file () {
	in_PLINK_sumstats=$( readlink -f $1 )
	ldsc_task_name=$2
	in_dir=$( dirname ${in_PLINK_sumstats} )	
	sub_dir_name=$( basename ${in_dir} )
	data_root=$( dirname $( dirname ${in_dir} ) )
	out_basename=$( basename ${in_PLINK_sumstats}  | cut -d '.' -f1-3 )
	echo "${data_root}/ldsc/${ldsc_task_name}/${sub_dir_name}/${out_basename}"
}
repo_dir=$(dirname $(dirname $(readlink -f $0)))


if [ $# -lt 1 ] ; then show_usage_and_exit ; fi
in_PLINK_sumstats=$( readlink -f $1 )

if [ $# -gt 1 ] ; then
	out_file_munged=$2
else
	out_file_munged=$( configure_out_file ${in_PLINK_sumstats} "munged" )
fi

if [ $# -gt 2 ] ; then
	out_file_h2=$3
else
	out_file_h2=$( configure_out_file ${in_PLINK_sumstats} "h2" )
fi

if [ $# -gt 3 ] ; then
	in_type=$4
elif [ $( echo $( basename ${in_PLINK_sumstats} ) | rev | cut -c -16 | rev ) == ".logistic.hybrid" ] ; then
	in_type="logistic"
elif [ $( echo $( basename ${in_PLINK_sumstats} ) | rev | cut -c -7 | rev )  == ".linear" ] ; then
	in_type="linear"
else 
	echo "Failed to infer the input type. Please pass in_type argument" >&2
	show_usage_and_exit
fi

in_name=$( basename ${in_PLINK_sumstats} )

if [ ! -f ${out_file_munged}.log ] ; then
	if [ ! -d $(dirname ${out_file_munged}) ] ; then mkdir -p $(dirname ${out_file_munged}) ; fi
	bash ${repo_dir}/helpers/ukb_ldsc_munge_wrapper.sh ${out_file_munged} ${in_name} ${in_type} ${in_PLINK_sumstats}
fi
echo ${out_file_munged}.log

if [ ! -f ${out_file_h2}.log ] ; then
	if [ ! -d $(dirname ${out_file_h2}) ] ; then mkdir -p $(dirname ${out_file_h2}) ; fi
	bash ${repo_dir}/helpers/ukb_ldsc_h2_helper.sh ${out_file_h2} ${out_file_munged}
fi
echo ${out_file_h2}.log



