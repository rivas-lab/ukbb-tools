#!/bin/bash
#SBATCH --job-name=imp_vqc
#SBATCH --output=logs/imp_vqc.%A.out
#SBATCH  --error=logs/imp_vqc.%A.err
#SBATCH --nodes=1
#SBATCH --cores=6
#SBATCH --mem=60000
#SBATCH --time=2:00:00
#SBATCH -p mrivas
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )


imp_dir="/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined"
pfile="${imp_dir}/pgen/ukb24983_ukb24983_cal_hla_cnv_imp"
vqc_prefix="${imp_dir}/variant_qc/ukb24983_ukb24983_cal_hla_cnv_imp"

plink_opts="--memory ${mem} --threads ${cores}"

plink2 ${plink_opts} --pfile ${pfile} vzs --out ${vqc_prefix} \
    --hardy zs midp cols=chrom,pos,ref,alt,ax,gcounts,hetfreq,sexaf,p
mv ${vqc_prefix}.log ${vqc_prefix}.hardy.log

plink2 ${plink_opts} --pfile ${pfile} vzs --out ${vqc_prefix} \
    --freq zs cols=chrom,pos,ref,alt,altfreq,machr2,nobs
mv ${vqc_prefix}.log ${vqc_prefix}.freq.log

plink2 ${plink_opts} --pfile ${pfile} vzs --out ${vqc_prefix} \
    --geno-counts zs cols=chrom,pos,ref,alt,homref,refalt,homalt1,altxy,hapref,hapalt,missing,nobs
mv ${vqc_prefix}.log ${vqc_prefix}.geno-counts.log
