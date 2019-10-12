#!/bin/bash
#SBATCH --job-name=maf1_no_cal
#SBATCH --output=logs/maf1_no_cal.%A.out
#SBATCH  --error=logs/maf1_no_cal.%A.err
#SBATCH --nodes=1
#SBATCH --cores=6
#SBATCH --mem=60000
#SBATCH --time=0:30:00
#SBATCH -p mrivas
set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

make_maf1_no_cal () {
    # c: chromosome
    #    it should be in [1-22],X,XY

    local c=$1
    local mem=$2
    local cores=$3

    local plink_opts="--memory ${mem} --threads ${cores}"
    local dir="/oak/stanford/groups/mrivas/ukbb/24983/imp/pgen"
    local out="${dir}/maf1_no_cal/ukb24983_imp_chr${c}_v3_maf1_no_cal"
    
    # exclude the variants on array and generate PLINK 2.0 files
    if [ ! -f ${out}.pgen.log ] ; then    
        cat ${dir}/maf1/ukb24983_cal_v2_hg19_imp_v3_maf1.join.tsv \
        | tail -n+2 | cut -f3 \
        | plink2 ${plink_opts} \
            --pfile ${dir}/maf1/ukb24983_imp_chr${c}_v3_maf1 vzs \
            --exclude /dev/stdin \
            --out ${out} \
            --make-pgen vzs        

        mv ${out}.log ${out}.pgen.log    
    fi 

    # prepare PLINK 1.9 files
    if [ ! -f ${out}.bed.log ] ; then    
        plink2 ${plink_opts} \
            --pfile ${out} vzs \
            --out ${out} \
            --make-bed

        mv ${out}.log ${out}.bed.log
    fi
}

chr=$1
make_maf1_no_cal ${chr} ${mem} ${cores}
