#!/bin/bash
set -beEuo pipefail

sqc_d="/oak/stanford/groups/mrivas/ukbb24983/sqc"
sqc_v="20200828"
pca_v="pca1_refinement"

pop=$1
shift 1

bash ../pca/1_pca.sh --nCores 6 --memory 40000 \
$([ "${pop}" == "e_asian" ] && echo "--EIG" || echo "") \
--keep ${sqc_d}/population_stratification_w24983_${sqc_v}/${pca_v}/ukb24983_${pop}.phe \
${sqc_d}/population_stratification_w24983_${sqc_v}/${pca_v}/ukb24983_${pop} \
$@

exit 0

##########################################
# job submission commands
##########################################

ml load resbatch

for pop in white_british ; do sbatch -p mrivas --qos=high_p --nodes=1 --mem=300000 --cores=10 --time=7-0:00:00 --job-name=pca.${pop} --output=logs/pca.${pop}.%A.out --error=logs/pca.${pop}.%A.err 2_pca1_refinement.sh ${pop} --nCores 10 --memory 300000 ; done

for pop in 'non_british_white' 'african' 's_asian' 'e_asian' ; do sbatch -p mrivas --qos=high_p --nodes=1 --mem=100000 --cores=6 --time=7-0:00:00 --job-name=pca.${pop} --output=logs/pca.${pop}.%A.out --error=logs/pca.${pop}.%A.err 2_pca1_refinement.sh ${pop} --nCores 6 --memory 100000 ; done

for pop in 'e_asian' ; do sbatch -p mrivas --qos=high_p --nodes=1 --mem=60000 --cores=6 --time=7-0:00:00 --job-name=pca.${pop} --output=logs/pca.${pop}.%A.out --error=logs/pca.${pop}.%A.err 2_pca1_refinement.sh ${pop} --nCores 6 --memory 60000 ; done
