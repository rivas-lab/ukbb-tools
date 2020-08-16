#!/bin/bash
set -beEuo pipefail

pop=$1

threads=6
memory="60G"

prefix='ukb24983_v2_hg19'
variant_type='array-combined'
freeze_v='20200815'
p_thr='1e-3'

data_d="/oak/stanford/groups/mrivas/ukbb24983/${variant_type}/gwas"
out_tsv="${data_d}/freeze/${freeze_v}/${prefix}.${pop}.${variant_type}.glm.${freeze_v}.p${p_thr}.tsv"

filtered_d="/scratch/groups/mrivas/ukbb24983/array-combined/gwas/freeze/${pop}"

    
echo '#CHROM POS Variant_ID GBE_ID population REF ALT A1 OBS_CT BETA SE P' | tr ' ' '\t' > ${out_tsv}

find ${filtered_d} -type f -size +68c | parallel "cat {} | grep -v '^#'" \
    | sort -S ${memory} --parallel=${threads} -k1,1V -k2,2n -k3,3V -k4,4V >> ${out_tsv}

bgzip --compress-level 9 -@ ${threads} -f ${out_tsv}
tabix -s 1 -b 2 -e 2 ${out_tsv}.gz

echo ${out_tsv}.gz

exit 0
# job submission command
for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' 'metal' ; do echo $pop ; sbatch -p mrivas --qos=high_p --nodes=1 --mem=64000 --cores=6 --time=2-0:00:00 --job-name=freeze_2b_combine.${pop} --output=logs/freeze_2b_combine.${pop}.%A.out --error=logs/freeze_2b_combine.${pop}.%A.err --dependency=afterok:5920082 freeze_2b_combine.sh ${pop} ; done
