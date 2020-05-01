#!/bin/bash
set -beEuo pipefail

_default_ldsc_path="$OAK/users/$USER/repos/bulik/ldsc"
_default_ld_scores="$OAK/projects/h2-estimation/private_output/ukb_ld_scores/TWINSUK/"
_default_tmp_dir="$TMPDIR"

# parse args
show_usage_and_exit () {
	echo "usage: $0 out munged1 munged2 [ld_scores ($_default_ld_scores)] [ldsc_path ($_default_ldsc_path)]" >&1 
	exit 1 
}
if [ $# -lt 3 ] ; then show_usage_and_exit ; fi
out=$1; munged1=$(readlink -f $2); munged2=$(readlink -f $3)
if [ $# -gt 3 ] ; then ld_scores=$4; else ld_scores=$_default_ld_scores; fi
if [ $# -gt 4 ] ; then ldsc_path=$5; else ldsc_path=$_default_ldsc_path; fi

# tmp_dir
tmp_dir=$(mktemp -d -p $_default_tmp_dir); echo "tmp_dir is $tmp_dir"
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

# helpfer functions
ldsc_ldsc () {
	ldsc_munged1=${1%.sumstats.gz}.sumstats.gz ;
	ldsc_munged2=${2%.sumstats.gz}.sumstats.gz ;
	ldsc_out=$3 ;

	# wrapper for ldsc.py 

	python $ldsc_path/ldsc.py \
	--rg ${ldsc_munged1},${ldsc_munged2} \
	--ref-ld-chr $ld_scores \
	--w-ld-chr $ld_scores \
	--out $ldsc_out 
}

if [ ! -d "$(dirname ${out})" ] ; then mkdir -p $(dirname ${out}) ; fi

out_basename=$(basename ${out} .log)

# run LDSC scripts
ldsc_ldsc ${munged1} ${munged2} ${tmp_dir}/${out_basename} 

# write the results
cp $tmp_dir/${out_basename}.log $(dirname ${out})/${out_basename}.log

echo "Results are written to:  ${out}.log"
