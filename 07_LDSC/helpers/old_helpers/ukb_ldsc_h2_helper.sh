#!/bin/bash
set -beEuo pipefail

_default_ldsc_path="$OAK/users/$USER/repos/bulik/ldsc"
_default_ld_scores="$OAK/projects/h2-estimation/private_output/ukb_ld_scores/TWINSUK/"
_default_tmp_dir="$TMPDIR"

# parse args
show_usage_and_exit () {
	echo "usage: $0 out munged [ld_scores ($_default_ld_scores)] [ldsc_path ($_default_ldsc_path)]" >&1 
	exit 1 
}
if [ $# -lt 2 ] ; then show_usage_and_exit ; fi
out=$1; munged=$(readlink -f $2); 
if [ $# -gt 2 ] ; then ld_scores=$3; else ld_scores=$_default_ld_scores; fi
if [ $# -gt 3 ] ; then ldsc_path=$4; else ldsc_path=$_default_ldsc_path; fi

# tmp_dir
tmp_dir=$(mktemp -d -p $_default_tmp_dir); echo "tmp_dir is $tmp_dir"
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

# helpfer functions
ldsc_ldsc () {
	ldsc_munged=${1%.sumstats.gz}.sumstats.gz ;
	ldsc_out=$2 ;

	# wrapper for ldsc.py 

	python $ldsc_path/ldsc.py \
	--h2 ${ldsc_munged} \
	--ref-ld-chr $ld_scores \
	--w-ld-chr $ld_scores \
	--out $ldsc_out 
}

out_basename=$(basename ${out} .log)

# run LDSC scripts
ldsc_ldsc ${munged} ${tmp_dir}/${out_basename} 

# write the results
cp $tmp_dir/${out_basename}.log $(dirname ${out})/${out_basename}.log

