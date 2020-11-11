#!/bin/bash
set -beEuo pipefail

pop=$1

prefix='ukb24983_v2_hg19'
variant_type='array-combined'
freeze_v='20201106'

data_d="/oak/stanford/groups/mrivas/ukbb24983/${variant_type}/gwas"
out_tar="${data_d}/freeze/${freeze_v}/${prefix}.${pop}.${variant_type}.glm.${freeze_v}.tar"

if [ ! -d ${data_d}/freeze/${freeze_v} ] ; then mkdir -p ${data_d}/freeze/${freeze_v}/ ; fi

find -L ${data_d}/current/${pop} -type f \
| while read f ; do readlink -f $f ; done \
| tar --transform 's/.*\///g' --group root --owner root -cvf ${out_tar} -T /dev/stdin

exit 0

# run in interactive bash session
for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' 'metal' ; do bash freeze_1_tar.sh ${pop} ; done

# submit as SLURM jobs
for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' 'metal' ; do sbatch -p mrivas --nodes=1 --mem=6000 --cores=1 --time=2-0:00:00 --job-name=freeze_tar.${pop} --output=logs/freeze_tar.${pop}.%A.out --error=logs/freeze_tar.${pop}.%A.err freeze_1_tar.sh ${pop} ; done
