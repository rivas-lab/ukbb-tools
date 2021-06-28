#!/bin/bash
#SBATCH --job-name=mrpmm_bin
#SBATCH --output=mrpmm_logs/mrpmm_bin.%A_%a.out
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=200000
#SBATCH --time=2-00:00:00

# define functions
usage () {
    echo "$0: MRPMM script for BIN phenotype exome data"
    echo "usage: sbatch -p <partition(s)> --array=1-<number of array jobs> $0 start_idx (inclusive) output_folder"
    echo "e.g. sbatch -p normal,owners --array=1-1000 $0 1 /path/to/output_folder"
}

# get core and memory settings from the header -- passed to gwas script below
cores=$(              cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(                cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )
log_dir=$( dirname $( cat $0 | egrep '^#SBATCH --output=' | awk -v FS='=' '{print $NF}' ))

# check number of command line args and dump usage if that's not right
if [ $# -ne 1 ] ; then usage >&1 ; exit 1 ; fi

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load python/3.6.1

_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname = $(hostname) SLURM_JOBID = ${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID = ${_SLURM_ARRAY_TASK_ID}" >&1

# get phenotypes to run
start_idx=$1
output_folder="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrpmm_exome"
this_idx=$_SLURM_ARRAY_TASK_ID

# read in the actual thing
GBE_ID=$(awk -F'\t' '{print $1}' mrp200kexomeresults.tsv | tail -n +2 | egrep -v QT_FC | egrep -v INI | sort -u | awk -v start_idx=$start_idx -v this_idx=$this_idx 'NR==(start_idx + this_idx - 1) {print $1}')
grep $GBE_ID mrp200kexomeresults.tsv | awk -F'\t' '{print $5}' > $output_folder/${GBE_ID}_genes.tsv
   
echo -e "$GBE_ID" >&1

rm $output_folder/$GBE_ID.tmp.txt
touch $output_folder/$GBE_ID.tmp.txt

POP="metal"
PATH_TO_FILE=$(find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/ -name "*.$GBE_ID.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep $POP | grep -v plink)
echo -e "$PATH_TO_FILE\t$POP\t$GBE_ID\tTRUE" >> $output_folder/$GBE_ID.tmp.txt
sed -i '1s/^/path\tstudy\tpheno\tR_phen\n/' $output_folder/$GBE_ID.tmp.txt

cat $output_folder/$GBE_ID.tmp.txt

for gene in $(cat $output_folder/${GBE_ID}_genes.tsv); do
    awk -F'\t' -v gene=$gene '{if ($4 == gene) {print $1}}' ../13_mrp/ref/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv > $output_folder/${GBE_ID}_${gene}_variants.tsv
    sed -i '1s/^/V\n/' $output_folder/${GBE_ID}_${gene}_variants.tsv
    python3 mrpmm.py --variants $output_folder/${GBE_ID}_${gene}_variants.tsv --phenotypes $output_folder/$GBE_ID.tmp.txt --metadata_path ../13_mrp/ref/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --out_folder $output_folder --C 1 2 3 4 5 --se_thresh 0.2 --maf_thresh 0.01 --variant_filters ptv pav
    rm $output_folder/${GBE_ID}_${gene}_variants.tsv
done

rm $output_folder/$GBE_ID.tmp.txt
rm $output_folder/${GBE_ID}_genes.tsv 
