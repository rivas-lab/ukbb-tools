#!/bin/bash
#SBATCH --job-name=mrp_rv_ma_array
#SBATCH --output=mrp_logs/mrp_rv_ma_array.%A_%a.out
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=32000
#SBATCH --time=04:00:00
#SBATCH -p normal,owners,mrivas

# define functions
usage () {
    echo "$0: MRP script to run rare-variant aggregation meta-analysis on array data"
    echo "usage: sbatch --array=1-<number of array jobs> $0 start_idx (inclusive) output_folder"
    echo "e.g. sbatch --array=1-1000 $0 1 /path/to/output_folder"
    echo '  You may check the status of the job (which jobs are finished) using the array-job module:'
    echo '  $ ml load array-job'
    echo '  $ array-job-find_ok.sh rerun_logs'
}

# get core and memory settings from the header -- passed to gwas script below
cores=$(              cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(                cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
log_dir=$( dirname $( cat $0 | egrep '^#SBATCH --output=' | awk -v FS='=' '{print $NF}' ))

# check number of command line args and dump usage if that's not right
if [ $# -ne 2 ] ; then usage >&1 ; exit 1 ; fi

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load python/2.7.13

# job start header (for use with array-job module)
_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&1

# get phenotypes to run
start_idx=$1
output_folder=$2
this_idx=$_SLURM_ARRAY_TASK_ID

min_N_count=10
GBE_ID=$(cat ../05_gbe/phenotype_info.tsv | awk -v min_N=${min_N_count} 'NR > 1 && $7 >= min_N' | egrep -v MED | awk -v start_idx=$start_idx -v this_idx=$this_idx 'NR==(start_idx + this_idx - 1) {print $1}' )

echo -e "path\tstudy\tpheno\tR_phen" > $output_folder/$GBE_ID.tmp.txt
echo -e "$GBE_ID" >&1

for POP in white_british african e_asian s_asian non_british_white; do
    lines=$(find /oak/stanford/groups/mrivas/ukbb24983/cal/gwas -name "*.$GBE_ID.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep $POP | wc -l)
    if [ $lines -eq 1 ]; then
        PATH_TO_FILE=$(find /oak/stanford/groups/mrivas/ukbb24983/cal/gwas -name "*.$GBE_ID.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep $POP)
        echo -e "$PATH_TO_FILE\t$POP\t$GBE_ID\tTRUE" >> $output_folder/$GBE_ID.tmp.txt
    fi
done

cat $output_folder/$GBE_ID.tmp.txt

/share/software/user/open/python/2.7.13/bin/python mrp_production.py --file $output_folder/$GBE_ID.tmp.txt --R_study independent similar --R_var independent similar --variants ptv pav --metadata_path /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep.tsv --out_folder $output_folder

rm $output_folder/$GBE_ID.tmp.txt
