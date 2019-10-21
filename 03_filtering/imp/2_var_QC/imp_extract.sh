#!/bin/bash
#SBATCH --job-name=imp
#SBATCH --output=logs/imp.%A.out
#SBATCH  --error=logs/imp.%A.err
#SBATCH --nodes=1
#SBATCH --cores=6
#SBATCH --mem=60000
#SBATCH --time=1-00:00:00
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

imp_extract () {
    # c: chromosome
    #    it should be in [1-22],X,XY

    local c=$1
    local mem=$2
    local cores=$3
    local tmp_dir=$4

    local plink_opts="--memory ${mem} --threads ${cores}"
    local imp_dir="/oak/stanford/groups/mrivas/ukbb/24983/imp"
    
    local var_lst="${imp_dir}/mfi/ukb_mfi_v3.info.maf.biallelic.noncal.tsv.zst"
    local tmp_var_lst=${tmp_dir}/extract.lst
    
    local in_prefix="${imp_dir}/pgen/ukb24983_imp_chr${c}_v3"
    local out_prefix="${imp_dir}/pgen/var_QC/ukb24983_imp_chr${c}_v3_QC"
    
    zstdcat ${var_lst} | egrep "^${c}:" | cut -f1 > ${tmp_var_lst}

    plink2 ${plink_opts} \
        --pfile   ${in_prefix} vzs \
        --extract ${tmp_var_lst} \
        --geno 0.01 \
        --hwe 1e-7 midp \
        --out     ${out_prefix} \
        --make-pgen vzs 

    mv ${out_prefix}.log ${out_prefix}.pgen.log 

    plink2 ${plink_opts} \
        --pfile ${out_prefix} vzs \
        --out ${out_prefix} \
        --make-bed

    mv ${out_prefix}.log ${out_prefix}.bed.log
}

chr=$1

imp_extract ${chr} ${mem} ${cores} ${tmp_dir}

