#!/bin/bash
#SBATCH --job-name=pgen
#SBATCH --output=logs/pgen.%A.out
#SBATCH  --error=logs/pgen.%A.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=40000
#SBATCH --time=2-0:00:00
# SBATCH -p mrivas,normal
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

out="/oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined/pgen_v2/ukb24983_ukb24983_cal_hla_cnv_imp"

plink_mem=$( perl -e "print(int(${mem} * 0.8))" )
plink_opts="--memory ${plink_mem} --threads ${cores}"
plink2 ${plink_opts} --bfile ${out} --keep-allele-order --make-pgen vzs --out ${out}

mv "${out}.log" "${out}.pgen.log"
