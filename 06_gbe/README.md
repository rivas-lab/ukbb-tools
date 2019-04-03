# Processing data for the Global Biobank Engine — a pipeline

This guide will explain the usage of the `gbe.sh` script in this folder. The purpose of this script is to conveniently process and conduct preliminary analyses (compatible with GBE) on a table of phenotype definitions as defined during a session (an example of which can be found [here](https://github.com/rivas-lab/ukbb-tools/blob/master/phenotyping/example_phenotyping_session.tsv)). This includes:

 - Processing phenotypes into PLINK-compatible specifications (`.phe` files)
 - Running GWAS on the set of directly genotyped variants in UK BioBank

In the future, this may be expanded to other analyses such as GWAS on imputed data, additional genotype models, GWAS under multiple populations, etc.

## Description and Implementation

This is a shell script, designed for use as a parallelized array job with `sbatch --array`. For those unfamiliar with the Slurm management system, [here](https://www.sherlock.stanford.edu/docs/user-guide/running-jobs/) is a description of how job submission works on Sherlock, and [here](https://slurm.schedmd.com/job_array.html) is some documentation on array jobs.

In short, the idea for this script is to iterate over each line in the input table, define the corresponding `.phe` file, then conduct analysis. 

In order to accomplish this, the script will do the following steps within each array instance:

 - Call [`tsv_to_phenos.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/phenotyping/scripts/tsv_to_phenos.py) using the `--only-this-row` option
 - Call [`gwas.py`](https://github.com/rivas-lab/ukbb-tools/blob/master/gwas/gwas.py) with the `--run-array` and `--run-now` options. The remaining flags are set to defaults used with GBE (namely, running analyses only on White British unrelated individuals).
 
_Note:_ The `--run-now` option will immediately run a GWAS (on directly genotyped data) directly on whichever node you get from the array job (as opposed to submitting a separate batch job, with parameters specified as arguments to `gwas.py`). This means you'll have to ensure the job specification for `gbe.sh` has enough compute (memory, time, etc.) to run this computation on each phenotype. The default parameters (16GB memory, 1 day of computation time) should suffice for most phenotypes, but in the event that the job fails, you will have to resubmit with more resources.

## Usage

Prior to running this script, you will need to specify an input table (line 15). The script will automatically go to this table and use appropriate columns via the table [here](https://docs.google.com/spreadsheets/d/1d4w4A8takvPxpHoUFXoNjj3a3QZLc-oHQiMaA5eElRg/edit?usp=sharing).

Then, you're ready to run your analysis. Run the following command where `path_to_tsv_file` is the string you wrote on line 15:

`sbatch --array=1-$(wc -l path_to_tsv_file) gbe.sh`
