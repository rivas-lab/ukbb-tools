#!/bin/bash
set -beEuo pipefail

# ml load R/3.6 gcc

# in_f=$1
in_f=/scratch/groups/mrivas/public_data/summary_stats/finngen_r3/summary_stats_hg19_plink/finngen_r3_AB1_ARTHROPOD.hg19.plink.gz

out_f=$(dirname $(dirname $in_f))/summary_stats_hg19_ldsc_munge/$(basename $in_f .plink.gz)

if [ ! -d $(dirname $out_f) ] ; then mkdir -p $(dirname $out_f) ; fi

if [ ! -f ${out_f}.log ] || [ ! -f ${out_f}.sumstats.gz ] ; then
    bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/helpers/ldsc_munge.sh --scratch --merge $in_f $out_f
fi

exit 0
###########################
# usage

find /scratch/groups/mrivas/public_data/summary_stats/finngen_r3/summary_stats_hg19_plink  -name "*.hg19.plink.gz"  | sort > 6_ldsc_munge.input.lst

sbatch -p mrivas,owners,normal --nodes=1 --mem=8000 --cores=1 --time=1:00:00 --job-name=ldsc_munge --output=logs/ldsc_munge.%A_%a.out --error=logs/ldsc_munge.%A_%a.err --array=1-901 /oak/stanford/groups/mrivas/users/ytanigaw/repos/yk-tanigawa/resbatch/parallel-sbatch.sh 6_ldsc_munge.sh 6_ldsc_munge.input.lst 2
