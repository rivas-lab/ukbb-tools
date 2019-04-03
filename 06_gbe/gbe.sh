#!/bin/bash
 
#SBATCH  --job-name=phenos
#SBATCH    --output=gbe.%A_%a.out
#SBATCH       --mem=16000
#SBATCH      --time=1-00:00:00
#SBATCH --partition=normal,owners

ml load plink2

# step 0: identify phenotype for processing
pheno_index=$(expr ${SLURM_ARRAY_TASK_ID} - 1)

# step 1: process phenotypes from input table
tsv_in="../02_phenotyping/tables/ukb_20190327.tsv"
gbe_input_tsv="../02_phenotyping/tables/gbe_sh_input_params.tsv"
gwasOutDir="/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/"

# look for table in reference above, throw error if it isn't
tsv_name="$(basename $tsv_in)"
relevant_row="$(grep $tsv_name $gbe_input_tsv)"
if [ -z $relevant_row ]; then
    echo "Could not find input table ${tsv_name} in reference table ${gbe_input_tsv}!"
    exit 2
fi

# automatically finds **zero-indexed column ids** for the below from gbe_input_tsv:
nameCol="$(cut -d' ' -f2 <<< $relevant_row)" # GBE ID
fieldCol="$(cut -d' ' -f3 <<< $relevant_row)" # Source UK Biobank Field ID (e.g. 21001, body mass index)
tableCol="$(cut -d' ' -f4 <<< $relevant_row)" # Source UK Biobank Table ID (e.g. 9797)
caseCol="$(cut -d' ' -f5 <<< $relevant_row)"  # Binary case codes
ctrlCol="$(cut -d' ' -f6 <<< $relevant_row)"  # Binary control codes
exclCol="$(cut -d' ' -f7 <<< $relevant_row)"  # Quantitative values to mark as missing
orderCol="$(cut -d' ' -f8 <<< $relevant_row)" # Order of quantitative values (least to greatest) in categorical fields 
descCol="$(cut -d' ' -f9 <<< $relevant_row)"   # String description of input phenotype (e.g. "Standing_height")

tableID="$(awk -F'\t' -v row=$SLURM_ARRAY_TASK_ID -v col=$tableCol -v h=1 '(NR==(row+h)){print $(col+1)}' $tsv_in)"
gwasOut="$(find $gwasOutDir -type d -name $tableID)"

# TODO: account for structure in phenotypedata directory due to basket id 
#      (this will likely have to be passed as a new argument)

# here's the command for phenotyping
python ../02_phenotyping/scripts/tsv_to_phenos.py  --tsv $tsv_in --only-this-row $pheno_index \
                                                   --name $nameCol --desc $descCol \
                                                   --field $fieldCol --table $tableCol \
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

# WARNING: this might colossally fail if we have duplicate pheno filenames

# find the phenotype file we just made
pheDir="/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata"
gbeId="$(awk -F'\t' -v row=$SLURM_ARRAY_TASK_ID -v col=$nameCol -v h=1 '(NR==(row+h)){print $(col+1)}' $tsv_in )"
echo $gbeId
pheFile=$( ls ${pheDir}/*/${tableID}/${gbeId}.phe )

# prep
mkdir -p ${gwasOutDir}/logs
logDir=`pwd`

# run gwas
python ../04_gwas/gwas.py --run-array --run-now --pheno $pheFile --out $gwasOut \
                          --population white_british --log-dir $logDir

# move log file
for type in genotyped exome; do 
    if [ -f ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ]; then
        mv ${gwasOutDir}/ukb24983_v2.${gbeId}.${type}.log ${gwasOutDir}/logs
    fi
done

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
