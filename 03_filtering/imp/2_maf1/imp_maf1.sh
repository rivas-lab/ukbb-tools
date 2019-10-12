#!/bin/bash
#SBATCH --job-name=imp_maf1
#SBATCH --output=logs/imp_maf1.%A.out
#SBATCH  --error=logs/imp_maf1.%A.err
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

make_imp_maf1 () {
    # c: chromosome
    #    it should be in [1-22],X,XY

    local c=$1
    local mem=$2
    local cores=$3
    local tmp_dir=$4

    local plink_opts="--memory ${mem} --threads ${cores}"
    local dir="/oak/stanford/groups/mrivas/ukbb/24983/imp/pgen"

    plink2 ${plink_opts} \
        --pfile ${dir}/ukb24983_imp_chr${c}_v3 vzs --maf 0.01 \
        --out   ${tmp_dir}/ukb24983_imp_chr${c}_v3_maf1 \
        --make-just-pvar zs 
    
    Rscript extract_monoallelic_vars.R \
    ${tmp_dir}/ukb24983_imp_chr${c}_v3_maf1.pvar.zst \
    ${tmp_dir}/extract.lst

    plink2 ${plink_opts} \
        --pfile ${dir}/ukb24983_imp_chr${c}_v3 vzs --maf 0.01 \
        --extract ${tmp_dir}/extract.lst \
        --out   ${dir}/maf1/ukb24983_imp_chr${c}_v3_maf1 \
        --make-pgen vzs 

    mv ${dir}/maf1/ukb24983_imp_chr${c}_v3_maf1.log ${dir}/maf1/ukb24983_imp_chr${c}_v3_maf1.pgen.log 
}

chr=$1

make_imp_maf1 ${chr} ${mem} ${cores} ${tmp_dir}

