#!/bin/bash

#SBATCH --job-name=exome
#SBATCH   --output=exome.%A.out
#SBATCH    --error=exome.%A.err
#SBATCH --time=2-0:00:00
#SBATCH --qos=normal
#SBATCH -p owners,normal,mrivas
#SBATCH --nodes=1
#SBATCH --cores=2
#SBATCH --mem=4000
#SBATCH --mail-type=END,FAIL
#################
set -beEuo pipefail

for fid in 23161 23171 23162 23172 ; do 
    outd=/scratch/groups/mrivas/ukbb/24983/phenotypedata/2003422/current/download/ukb28249.${fid}
    if [ ! -d $outd ] ; then mkdir -p $outd ; fi
    bash /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/08_bulk_DL/ukbfetch_bulk_wrapper.sh \
        /oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/2003422/current/download/ukb28249.${fid}.bulk \
        $outd \
        /oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/2003422/current/download/ukb28249.key \
        10
done

