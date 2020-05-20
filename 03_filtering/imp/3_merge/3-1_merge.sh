#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=logs/merge.%A.out
#SBATCH  --error=logs/merge.%A.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=60000
#SBATCH --time=7-0:00:00
#SBATCH -p mrivas
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
source $(dirname $(readlink -f $0))/3-0_merge_misc.sh
##############

merge_list_tsv=$1
out=$2

plink_mem=$( perl -e "print(int(${mem} * 0.8))" )
plink_opts="--memory ${plink_mem} --threads ${cores}"
tmp_merge_list=${tmp_dir}/merge.lst

# prepare temp files for merge (assign short IDs)
cat ${merge_list_tsv} | while read f p ; do 
    echo "processing $f" >&2    
    prep_files $f $p $tmp_dir $(dirname $(readlink -f $0))/assign_short_names_for_merge.R
    echo ${tmp_dir}/$(basename $f) | sed -e "s/.zst$//g" | sed -e "s/.pvar$//g" >> ${tmp_merge_list}
done

cat ${tmp_merge_list} >&2

plink ${plink_opts} --merge-list ${tmp_merge_list} --make-bed --keep-allele-order --out ${out}

mv "${out}.bim" "${out}.shortnames.bim"
