#!/bin/bash
#SBATCH --job-name=bulk-dl
#SBATCH --output=bulk-dl.%A.out
#SBATCH  --error=bulk-dl.%A.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=30000
#SBATCH --time=7-00:00:00
#SBATCH -p mrivas

set -beEuo pipefail

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

fid=$1

ml load ukbb-showcase-utils

out_d="/scratch/groups/mrivas/ukbb24983/phenotypedata/2005693/41413/bulk/ukb2005693.41413.${fid}"
src="/oak/stanford/groups/mrivas/users/${USER}/repos/rivas-lab/ukbb-tools/08_bulk_DL/ukbfetch_bulk_wrapper.sh"

echo bash ${src} ${out_d}/$(basename ${out_d}).bulk ${out_d} ${out_d}/.ukbkey ${cores}
bash ${src} ${out_d}/$(basename ${out_d}).bulk ${out_d} ${out_d}/.ukbkey ${cores}

