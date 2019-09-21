#!/bin/bash
#SBATCH --job-name=split
#SBATCH --output=logs/split.%A.out
#SBATCH  --error=logs/split.%A.err
#SBATCH --nodes=1
#SBATCH --cores=6
#SBATCH --mem=60000
#SBATCH --time=1-0:00:00
#SBATCH -p mrivas
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

imp_dir="/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined"
pfile="${imp_dir}/pgen/ukb24983_ukb24983_cal_hla_cnv_imp"
out="${imp_dir}/split/ukb24983_ukb24983_cal_hla_cnv_imp"
split="/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20190809/split_array"

plink_opts="--memory ${mem} --threads ${cores}"

for s in val test train ; do
    plink2 ${plink_opts} --pfile ${pfile} vzs \
        --keep ${split}/ukb24983_white_british_${s}.phe \
        --make-pgen vzs --out ${out}.${s}
    mv ${out}.${s}.log ${out}.${s}.pgen.log
    
    plink2 ${plink_opts} --pfile ${out}.${s} vzs \
        --geno-counts zs cols=chrom,pos,ref,alt,homref,refalt,homalt1,altxy,hapref,hapalt,missing,nobs \
        --out ${out}.${s}
    mv ${out}.${s}.log ${out}.${s}.geno-counts.log

    plink2 ${plink_opts} --pfile ${out}.${s} vzs \
        --make-bed --out ${out}.${s}
    mv ${out}.${s}.log ${out}.${s}.bed.log
done
