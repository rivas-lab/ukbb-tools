#!/bin/bash
#SBATCH --job-name=LD_MAP
#SBATCH --output=logs/ldmap.%A.out
#SBATCH  --error=logs/ldmap.%A.err
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=64000
#SBATCH --time=2-00:00:00
#SBATCH -p mrivas,normal

SRCNAME=$(readlink -f $0)
SRCDIR=$(dirname ${SRCNAME})
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"

set -beEuo pipefail
############################################################

source /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/14_LD_map/14_LD_map_misc.sh
ml load plink plink2

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

#########################
memory=$(perl -e "print(int(${mem} * 0.8))")
# for plink, we allocate the 80% of the memory
bfile="/scratch/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv"
out_prefix="/scratch/groups/mrivas/ukbb24983/array_combined/ldmap_train/ukb24983_cal_hla_cnv.train"
keep_file="${tmp_dir}/train.phe"

# generate keep file from GWAS covar file
cat /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200313/ukb24983_GWAS_covar.20200313.phe \
| awk -v FS='\t' -v OFS='\t' '($4=="train"){print $1, $2}' > ${keep_file}

plink_opts="" # one can put maf threshold etc.

compute_ld_map_wrapper \
    "${bfile}" \
    "${out_prefix}" \
    "${cores}" "${memory}" \
    "${keep_file}" \
    "${plink_opts}"

# job finish footer (for use with array-job module)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2
