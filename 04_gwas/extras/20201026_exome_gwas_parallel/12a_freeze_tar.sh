#!/bin/bash
set -beEuo pipefail

pop=$1

variant_type='exome'
freeze_v='master_phe_20201002_exomeOQFE_20201110'

# data_d="/oak/stanford/groups/mrivas/ukbb24983/${variant_type}/gwas"
data_d="/scratch/groups/mrivas/ukbb24983/${variant_type}/gwas"
out_tar="${data_d}/freeze/${freeze_v}/${freeze_v}.${pop}.tar"

if [ ! -d ${data_d}/freeze/${freeze_v} ] ; then mkdir -p ${data_d}/freeze/${freeze_v}/ ; fi

find -L ${data_d}/master_phe_20201002_exomeOQFE/${pop} -type f \
| while read f ; do readlink -f $f ; done \
| tar --transform 's/.*\///g' --group root --owner root -cvf ${out_tar} -T /dev/stdin

exit 0
#################################

for pop in 'white_british' 'non_british_white' 'african' 's_asian' 'e_asian' 'related' 'others' 'metal' ; do
    sbatch -p mrivas --nodes=1 --mem=12000 --cores=2 --time=7-0:00:00 \
        --job-name=tar.${pop} --output=logs/tar.${pop}.%A.out --error=logs/tar.${pop}.%A.err \
        12a_freeze_tar.sh ${pop}
done

# Then upload to gdrive
for pop in e_asian s_asian african related others non_british_white ; do echo "[$(date +%Y%m%d-%H%M%S)] ${pop}" ; rclone copy /scratch/groups/mrivas/ukbb24983/exome/gwas/freeze/master_phe_20201002_exomeOQFE_20201110/master_phe_20201002_exomeOQFE_20201110.${pop}.tar gdrive://rivas-lab/ukbb24983/exome/freeze/master_phe_20201002_exomeOQFE_20201110 ; done

