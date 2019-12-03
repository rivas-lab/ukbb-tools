#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=logs/merge.%A.out
#SBATCH  --error=logs/merge.%A.err
#SBATCH --nodes=1
#SBATCH --cores=6
#SBATCH --mem=64000
#SBATCH --time=0-20:00:00
#SBATCH -p mrivas,normal
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

##############
tmp_dir_root=${LOCAL_SCRATCH}
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
##############

helper_src_dir="/oak/stanford/groups/mrivas/users/$USER/repos/rivas-lab/ukbb-tools/03_filtering/imp/3_merge"

source "${helper_src_dir}/3-0_merge_misc.sh"

##############

merge_list_tsv=$1
out="/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined_no_cnv/pgen/ukb24983_ukb24983_cal_hla_imp"

plink_opts="--memory 40000 --threads ${cores}"
tmp_merge_list=${tmp_dir}/merge.lst

# prepare temp files for merge (assign short IDs)
cat ${merge_list_tsv} | while read f p ; do 
    echo "processing $f" >&2
    prep_files $f $p $tmp_dir ${helper_src_dir}/assign_short_names_for_merge.R
    echo ${tmp_dir}/$(basename $f) | sed -e "s/.zst$//g" | sed -e "s/.pvar$//g" >> ${tmp_merge_list}
done

cat ${tmp_merge_list} >&2

plink ${plink_opts} --merge-list ${tmp_merge_list} --make-bed --keep-allele-order --out ${out}

mv "${out}.bim" "${out}.shortnames.bim"
