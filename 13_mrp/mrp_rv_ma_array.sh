#!/bin/bash
#SBATCH --job-name=mrp_rv_ma_array
#SBATCH --output=mrp_logs/mrp_rv_ma_array.%A_%a.out
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=32000
#SBATCH --time=04:00:00

# define functions
usage () {
    echo "$0: MRP script to run rare-variant aggregation meta-analysis on array data"
    echo "usage: sbatch -p <partition(s)> --array=1-<number of array jobs> $0 start_idx (inclusive) output_folder"
    echo "e.g. sbatch -p normal,owners --array=1-1000 $0 1 /path/to/output_folder"
}

# get core and memory settings from the header -- passed to gwas script below
cores=$(              cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(                cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
log_dir=$( dirname $( cat $0 | egrep '^#SBATCH --output=' | awk -v FS='=' '{print $NF}' ))

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

echo -e "path\tstudy\tpheno\tR_phen" > $output_folder/$GBE_ID.tmp.txt
echo -e "$GBE_ID" >&1

line_phenotype="$(awk -F'\t' -v GBE=${GBE_ID} '{if ($1 == GBE) {print}}' ../05_gbe/array-combined/phenotype_info.tsv)"

echo -e "white_british\nafrican\ne_asian\ns_asian\nnon_british_white\nrelated\nothers" > $GBE_ID.POP

N_GBE="$(echo $line_phenotype | cut -d' ' -f8)"
N_NBW="$(echo $line_phenotype | cut -d' ' -f9)"
N_AFR="$(echo $line_phenotype | cut -d' ' -f10)"
N_EAS="$(echo $line_phenotype | cut -d' ' -f11)"
N_SAS="$(echo $line_phenotype | cut -d' ' -f12)"
N_REL="$(echo $line_phenotype | cut -d' ' -f13)"
N_OTH="$(echo $line_phenotype | cut -d' ' -f14)"

echo -e "$N_GBE\n$N_AFR\n$N_EAS\n$N_SAS\n$N_NBW\n$N_REL\n$N_OTH" > $GBE_ID.POP_NUM

paste $GBE_ID.POP $GBE_ID.POP_NUM | while read POP NUM; do
    lines=$(find /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/ -name "*.$GBE_ID.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep $POP | wc -l)
    if [ $lines -eq 1 ] && [ $NUM -ge $min_N_count ]; then
        PATH_TO_FILE=$(find /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas/current/ -name "*.$GBE_ID.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep $POP)
        echo -e "$PATH_TO_FILE\t$POP\t$GBE_ID\tTRUE" >> $output_folder/$GBE_ID.tmp.txt
    fi
done

cat $output_folder/$GBE_ID.tmp.txt

/share/software/user/open/python/3.6.1/bin/python3 mrp_production.py --file $output_folder/$GBE_ID.tmp.txt --R_study independent similar --R_var independent similar --variants ptv pav --metadata_path /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --out_folder $output_folder

rm $output_folder/$GBE_ID.tmp.txt
rm $GBE_ID.POP_NUM
rm $GBE_ID.POP
