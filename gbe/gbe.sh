#!/bin/bash
 
#SBATCH  --job-name=phe2gbe
#SBATCH    --output=logs/gbe_pipeline.test.%A_%a.out
#SBATCH       --mem=24000
#SBATCH      --time=1-00:00:00
#SBATCH --partition=normal,owners

ml load plink2

# step 0: identify phenotype for processing
pheno_index=$(expr ${SLURM_ARRAY_TASK_ID} - 1)


# step 1: process phenotypes from input table
tsv_in="/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes/tables/on_github/ukb_20170818.tsv"

# provide **zero-indexed column ids** for the below:
nameCol=3 # GBE ID
fieldCol=6 # Source UK Biobank Field ID (e.g. 21001, body mass index)
tableCol=4 # Source UK Biobank Table ID (e.g. 9797)
caseCol=13  # Binary case codes
ctrlCol=14  # Binary control codes
exclCol=11  # Quantitative values to mark as missing
orderCol=12 # Order of quantitative values (least to greatest) in categorical fields 

# TODO: account for structure in phenotypedata directory due to basket id 
#      (this will likely have to be passed as a new argument)

# here's the command
python ../phenotyping/scripts/tsv_to_phenos.py  --tsv $tsv_in --only-this-row $pheno_index \
                                                --name $nameCol --field $fieldCol --table $tableCol \
                                                --case $caseCol --control $ctrlCol \
                                                --missing $exclCol --order $orderCol 
# extra options: --no-header, if the input TSV doesn't have a header row
#                --missing-is-control, to expand the control set to the entire population

# and the full readme for the script:
<<"COMMENT"
"""
  --tsv INPUT           input table from phenotyping session
  --no-header           flag if input tsv has no header
  --name NAME           column in tsv corresponding to phe file name (GBE ID)
  --field FIELD         column in tsv corresponding to UK Biobank Field ID
  --table TABLE         column in tsv corresponding to UK Biobank Table ID
  --case CASE           column in tsv corresponding to values for binary case
                        definitions
  --control CONTROL     column in tsv corresponding to values for binary
                        control definitions
  --missing-is-control  flag if missing values for binary traits should be
                        defined as controls
  --missing EXCLUDE     column in tsv corresponding to QT values considered as
                        missing data
  --order ORDER         column in tsv corresponding to order of values (least
                        to greatest) for QTs from categorical fields
  --only-this-row ONLYONE
                        (optional) flag to run only one (zero-indexed, not
                        including the header) row the input tsv.
"""
COMMENT



# step 2: run gwas on resulting phenotypes

# this finds the path to the newly defined phenotype (it's defined in make_phe.py)
# h is a hack -- needs to be 1/0 if input tsv has/lacks a header
# the plus ones are because awk is 1-indexed, and the provided columns are zero-indexed

# WARNING: this will colossally fail if we have duplicate pheno filenames, so let's not do that

# pheDir="/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata"
pheDir="/oak/stanford/groups/mrivas/dev-ukbb-tools/phenotypes"
gbeId="$(awk -F'\t' -v row=$SLURM_ARRAY_TASK_ID -v col=$nameCol -v h=1 '(NR==(row+h)){print $(col+1)}' $tsv_in )"
echo $gbeId
pheFile=$( find ${pheDir} -type f -name "${gbeId}.phe" )

# here are some more (programmable, i guess) parameters for the gwas:
gwasOutDir="/oak/stanford/groups/mrivas/dev-ukbb-tools/gwas/ukb_20170818"
mkdir -p $gwasOutDir
logDir=`pwd`

python ../gwas/gwas.py --run-array --run-now --pheno $pheFile --out $gwasOutDir \
                       --population white_british --log-dir $logDir

# and here's the readme for the gwas script:
<<"COMMENT"
"""
  --plink1              Flag to run GWAS with PLINK v1.9 instead of v2.0
  --run-array           Run GWAS on directly genotyped (array) data
  --run-imputed         Run GWAS on imputed data
  --pheno [PHENO [PHENO ...]]
                        Path to phenotype file(s)
  --out OUTDIR          Path to desired output *directory*. Summary stats will
                        be output according to phenotype name (derived from
                        passed file) and Rivas Lab specification for GBE.
                        Defaults to current working directory.
  --population POP      Flag to indicate which ethnic group to use for GWAS.
                        Must be one of all, white_british, e_asian, s_asian,
                        african
  --keep-related        Flag to keep related individuals in GWAS. Default is
                        to remove them.
  --batch-memory SB_MEM
                        For underlying batch job submission: Amount of memory
                        (in MB) to request. Default is 16000.
  --batch-time SB_TIME  For underlying batch job submission: Amount of time
                        (DD-HH:MM:SS) to request. Default is 24 hours.
  --batch-partitions [SB_PARTI [SB_PARTI ...]]
                        For underlying batch job submission: Compute partition
                        to submit jobs. Default is normal,owners.
  --log-dir LOG         Directory in which to place log files from this script
                        and its related SLURM jobs. Default is current working
                        directory.
  --run-now             Flag to run GWAS immediately, on the local machine,
                        rather than submitting a script with sbatch. Not
                        available for use with --run-imputed.
"""
COMMENT

echo "Completed."

# next, once the phenotypes are all defined, we'll have to update the reference documents
# (phenotype reference, master phenotype file, and icdinfo.txt)
