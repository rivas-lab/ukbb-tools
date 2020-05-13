#!/bin/bash
set -beEuo pipefail

_default_ldsc_path="$OAK/users/$USER/repos/bulik/ldsc"
_default_ld_scores="$OAK/projects/h2-estimation/private_output/ukb_ld_scores/TWINSUK/"
_default_tmp_dir="$TMPDIR"
_ldsc_input_script="$(dirname $(readlink -f $0))/make_ldsc_input_file_v2.R"

# parse args
show_usage_and_exit () {
	echo "$0 out_file in_name [in_type in_sumstats] [ld_scores ($_default_ld_scores)] [ldsc_path ($_default_ldsc_path)]" >&1 
	exit 1 
}
if [ $# -lt 2 ] ; then show_usage_and_exit ; fi
out_file=$1; in_name=$2;
if [ $# -gt 2 ] ; then
	if [ $# -lt 4 ] ; then show_usage_and_exit ; fi
	in_type=$3; in_sumstats=$4;
else
	in_type=""; in_sumstats="";
fi
if [ $# -gt 4 ] ; then ld_scores=$5; else ld_scores=$_default_ld_scores; fi
if [ $# -gt 5 ] ; then ldsc_path=$6; else ldsc_path=$_default_ldsc_path; fi

# tmp_dir
tmp_dir=$(mktemp -d -p $_default_tmp_dir); echo "tmp_dir is $tmp_dir"
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

# helpfer functions
ldsc_munge ()  {
	munge_sumstats=$1
	munge_out=$2
	ldsc_path=$3

	# wrapper for munge_sumstats.py
	python ${ldsc_path}/munge_sumstats.py \
	--sumstats ${munge_sumstats} \
	--N-col OBS_CT --a1 A1 --a2 A2 --snp ID \
	--signed-sumstats BETA,0 \
	--out ${munge_out}
}

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

if [ ! -d "$(dirname ${out_file})" ] ; then mkdir -p $(dirname ${out_file}) ; fi

ml load snpnet_yt
echo Rscript ${_ldsc_input_script} $tmp_dir/${in_name}_ldsc.tsv ${in_sumstats} ${ld_scores}
Rscript ${_ldsc_input_script} $tmp_dir/${in_name}_ldsc.tsv ${in_sumstats} ${ld_scores}

# run LDSC scripts
ldsc_munge $tmp_dir/${in_name}_ldsc.tsv $tmp_dir/${in_name}.munge ${ldsc_path}

# write the results
cp $tmp_dir/${in_name}.munge.log ${out_file%.gz}.log
cp $tmp_dir/${in_name}.munge.sumstats.gz "${out_file%.sumstats.gz}.sumstats.gz"

echo "Results are written to: ${out_file%.sumstats.gz}.sumstats.gz"