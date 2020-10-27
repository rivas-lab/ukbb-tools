#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --output=logs/merge.%A.out
#SBATCH  --error=logs/merge.%A.err
#SBATCH --nodes=1
#SBATCH --cores=12
#SBATCH --mem=60000
#SBATCH --time=7-0:00:00
#SBATCH -p mrivas
#SBATCH --qos=high_p
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

##############
# file locations
#data_d="/scratch/groups/mrivas/ukbb24983/exome-oqfe2020"
data_d="/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020"
out=${data_d}/ukb24983_exomeOQFE
##############

ml load plink plink2

##############
tmp_dir_root=${LOCAL_SCRATCH:=/tmp}
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT
##############

plink_mem=$( perl -e "print(int(${mem} * 0.8))" )
plink_opts="--memory ${plink_mem} --threads ${cores}"

tmp_merge_list=${tmp_dir}/merge.lst

for c in $(seq 1 22) X Y ; do
    bed=${data_d}/download/ukb23155_c${c}_b0_v1.bed
    fam=${data_d}/download/ukb23155_c${c}_b0_v1_s200632.fam
    bim=${data_d}/download/UKBexomeOQFE_chr${c}.bim

    echo ${bed} ${bim} ${fam}
done > ${tmp_merge_list}

cat ${tmp_merge_list} >&2

echo plink ${plink_opts} --merge-list ${tmp_merge_list} --make-bed --keep-allele-order --mac 1 --out ${out}
plink --silent ${plink_opts} --merge-list ${tmp_merge_list} --make-bed --keep-allele-order --mac 1 --out ${out}

mv "${out}.log" "${out}.bed.log"

echo plink2 ${plink_opts} --bfile ${out} --keep-allele-order --make-pgen vzs --out ${out}
plink2 --silent ${plink_opts} --bfile ${out} --keep-allele-order --make-pgen vzs --out ${out}

mv "${out}.log" "${out}.pgen.log"

