#!/bin/bash
set -beEuo pipefail

GBE_ID=$1

source 0_functions.sh

# constants

data_d="/scratch/groups/mrivas/ukbb24983/exome/gwas/master_phe_20201002_exomeOQFE"
pops=('white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others')
metal_src=/oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/18_metal/run_metal.sh

# main    
out_prefix="${data_d}/metal/ukb24983_exomeOQFE.${GBE_ID}"
if [ ! -f ${out_prefix}.metal.info.txt ] && [ ! -f ${out_prefix}.metal.tsv.gz ] ; then
    for pop in ${pops[@]} ; do
            
        echo ${data_d}/${pop}/ukb24983_exomeOQFE.${GBE_ID}.$(get_plink_suffix ${GBE_ID}).gz
    done | while read f ; do
        if [ -f $f ] ; then echo $f ; fi
    done | bash ${metal_src} --nCores 4 --assembly hg38 -o ${out_prefix} -f /dev/stdin
fi

exit 0
#############################
# instruction for job-submission on Sherlock

bash 6a_metal.sh INI50 # 20 min-ish

metal_trait_lst=6b_metal.trait.lst.gen.20201101-204502.lst

sbatch -p mrivas,normal --time=3:0:00 --mem=24000 --nodes=1 --cores=4 \
--job-name=metal --output=logs_scratch/metal.%A_%a.out --error=logs_scratch/metal.%A_%a.err \
--array=1-609 \
${parallel_sbatch_sh} 6a_metal.sh \
${metal_trait_lst} 3

