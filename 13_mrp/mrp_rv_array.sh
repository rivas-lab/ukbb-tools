#!/bin/bash
#SBATCH --job-name=mrp_rv_array
#SBATCH --output=mrp_logs/mrp_rv_array.%A_%a.out
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=16000
#SBATCH --time=01:00:00

# define functions
usage () {
    echo "$0: MRP script to run rare-variant aggregation for white british array data"
    echo "usage: sbatch -p <partition(s)> --array=1-<number of array jobs> $0 start_idx (inclusive) output_folder"
    echo "e.g. sbatch -p normal,owners --array=1-1000 $0 1 /path/to/output_folder"
}

# check number of command line args and dump usage if that's not right
if [ $# -ne 2 ] ; then usage >&1 ; exit 1 ; fi

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load python/3.6.1

_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&1

# get phenotypes to run
start_idx=$1
output_folder=$2
this_idx=$_SLURM_ARRAY_TASK_ID

min_N_count=100
GBE_ID=$(cat ../05_gbe/array-combined/phenotype_info.tsv | awk -v min_N=${min_N_count} 'NR > 1 && $8 >= min_N' | awk -v start_idx=$start_idx -v this_idx=$this_idx 'NR==(start_idx + this_idx - 1) {print $1}' )
POP="white_british"
echo $GBE_ID >&1
FILEPATH=$(find /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/ -name "*.$GBE_ID.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep $POP);

echo -e "path\tstudy\tpheno\tR_phen\n$FILEPATH\t$POP\t$GBE_ID\tTRUE" > $output_folder/$GBE_ID.tmp.txt;

cat $output_folder/$GBE_ID.tmp.txt

/share/software/user/open/python/3.6.1/bin/python3 mrp_production.py --file $output_folder/$GBE_ID.tmp.txt --R_var independent similar --variants ptv pav --metadata_path /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --out_folder $output_folder

rm $output_folder/$GBE_ID.tmp.txt
