#!/bin/bash
#SBATCH --job-name=conv2pgen
#SBATCH --output=logs/conv2pgen.%A.out
#SBATCH  --error=logs/conv2pgen.%A.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=60000
#SBATCH --time=2-0:00:00
#SBATCH -p mrivas
#SBATCH --qos=high_p

set -beEuo pipefail

ml load plink2

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

out="/scratch/groups/mrivas/ukbb24983/array_imp_combined_no_cnv/pgen_v2/ukb24983_cal_hla_imp"

plink_opts="--memory ${mem} --threads ${cores}"
plink2 ${plink_opts} --bfile ${out} --keep-allele-order --make-pgen vzs --out ${out}

mv "${out}.log" "${out}.pgen.log"
