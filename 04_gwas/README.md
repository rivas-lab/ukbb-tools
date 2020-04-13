# Genome-wide association studies with UK Biobank data

The "main" script in this directory, `gwas.py`, is designed to be a convenient wrapper for running genome-wide associations with PLINK using the lab's resources from the UK Biobank. This directory contains a plethora of options and double-wrapper scripts around `gwas.py`.

## Contents

1. [`04_gwas_misc.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/04_gwas_misc.sh): A file containing common functions used across the GWAS scripts. Loaded by running `source 04_gwas_misc.sh`.
2. [`gwas.py`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/gwas.py): The bread and butter of this directory. This script uses a Python command-line interface to create the recommended PLINK script and submits it to the cluster. For details, see [this section](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas#gwas-options).
3. [`gwas_freeze.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/gwas_freeze.sh): Freezes the current summary statistics, capping the results at a certain p-value threshsold in order to keep the file concise. Useful for GBE uploading.
- Inputs: 
  - `check_file`: Output from `check_gwas.sh` that aggregates the paths to the data in one place.
  - `prefix`: Prefix for GWAS summmary statistics outputs; e.g. `ukb24983_v2_hg19`.
  - `variant_type`: Usually, what follows `prefix` and `pop` in summary statistic output files; e.g. `genotyped`.
  - `freeze_v`: Freeze version. Suggested that this is input as `YYYYMMDD`.
  - `pop`: One of `white_british`, `non_british_white`, `african`, `s_asian`, or `e_asian`.  Default: `white_british`.
  - `p_val_threshold`: Float representing the maximum p-value that will be contained in the tarball. Default: `1e-3`.
  - `out_dir_root`: Where the tarball will go. Default: `/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/freeze/`.
- Outputs: A new tarball in `/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/freeze/[freeze_v]/[pop]`, and symlinks to this version within the folder.
- Usage: `gwas_freeze.sh check_file prefix variant_type freeze_v [pop] [p_val_threshold] [out_dir_root]`
- Verbose uage: `gwas_freeze.sh check_array_gwas.XXXXXXX.out ukb24983_v2_hg19 genotyped 20200412 e_asian 1e-5 /desired/path`
4. [`plink_updater.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/plink_updater.sh): Installs a version of PLINK2 on Sherlock given the version date.
- Inputs: `date`, a `YYYYMMDD`-formatted string dictating the date of the version that the script should pull down and install.
- Outputs: A new `plink2` version installed as both non-AVX2 and AVX2 versions on Sherlock. *NOTE*: In order to retain consistency between jobs, the new version is **NOT** made the default version immediately as a result of this script.
- Example usage: `bash plink_updater.sh 20200409`
5. [`run_array_hla_cnv.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/run_array_hla_cnv.sh): The default workhorse for running GWAS across phenotypes; as indicated by the name, it runs GWAS across the array, HLA, and CNV data.
- Inputs: None, explicitly. The array parameters (number of jobs and start index) are instead specified.
- Outputs: GWAS summary statistics for the `array-combined` (array, HLA, and CNV) dataset in the correct basket/table organized folders in `/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/`, and also symlinks in `/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/current`.
- Usage: `sbatch --array=1-[desired number of jobs, max 1000] run_array_hla_cnv.sh [start_idx] [pop]`
- Example usage: `sbatch --array=1-1000 run_array_hla_cnv.sh 1 african`
6. [`run_exome.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/run_exome.sh): The default workhorse for running GWAS across phenotypes; as indicated by the name, it runs GWAS across the exome data.
- Inputs: None, explicitly. The array parameters (number of jobs and start index) are instead specified.
- Outputs: GWAS summary statistics for the `exome-spb` (Regeneron calls) dataset in the correct basket/table organized folders in `/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/`, and also symlinks in `/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current`.
- Usage: `sbatch --array=1-[desired number of jobs, max 1000] run_exome.sh [start_idx] [pop]`
- Example usage: `sbatch --array=1-1000 run_exome.sh 1 african`
7. [`show_gwas.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/show_gwas.sh): Shows a summary-level GWAS statistics file over many individual phenotypes, filtering for those statistics that are below the p-value threshold specified.
- Inputs:
  - `check_file`: Output from [`check_gwas/check_gwas.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/check_gwas/check_gwas.sh) that aggregates the paths to the data in one place.
  - `p_val_threshold`: Float representing the maximum p-value that will be contained in the tarball. Default: `1e-4`.
- Outputs: A file that shows capped summary statistics for variants across phenotypes.
- Example usage: `bash show_gwas.sh check_array_gwas.XXXXXXX.out 1e-3 | less`
8. [`check_gwas`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/check_gwas): Contains scripts for checking the status of array and exome GWAS runs. The scripts parallelly generate a summary of `.phe` files, corresponding GWAS summary statistics, and a count of the number of total and non-NA lines in the summary statistics.
9. [`extras`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/extras): This directory contains various miscellaneous scripts and files describing our `extras` GWAS - from Biobank Japan, MVP, and the Broad collaboration.
10. [`rerun_logs`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/rerun_logs): This directory contains the log files from the SLURM submissions of the `run_gwas.sh`-type files. None of the logs will show up on the GitHub due to a nifty `.gitignore`.

## Should these be moved?

`flipfix_and_liftOver.py`
`flipfix`
`ref_alt` 
`test-flipfix_and_liftOver.sh`

## Pipelines and Workflows

### GWAS options

A full list of options can be obtained by running `python gwas.py -h`, and this output is replicated here for reference:

```
(base) user$ python gwas.py -h
usage: gwas.py [-h] [--plink1] [--run-array] [--run-exome] [--run-exome-gatk]
               [--run-imputed] [--run-cnv] [--run-cnv-burden] [--run-hla]
               [--run-array-imp-combined] [--run-array-combined] --pheno PHENO
               --out OUTDIR [--population POP] [--keep KEEP]
               [--remove-add RMADD] [--keep-related] [--cores CORES]
               [--memory MEM] [--batch-time SB_TIME]
               [--batch-partitions [SB_PARTI [SB_PARTI ...]]] [--log-dir LOG]
               [--run-now] [--sex-div SEX_DIV] [--keep-sex KEEP_SEX]
               [--keep-sex-file KEEP_SEX_FILE] [--include-x INCLUDE_X]
               [--additional-plink-opts [PLINK_OPTS [PLINK_OPTS ...]]]

A script for running GWAS with UK Biobank data (array/exome/imputed/cnv genotypes) using PLINK.

Author: Matthew Aguirre (SUNET: magu)

(Updated 10/14/2019 by E Flynn for sex-div analysis)
(Updated 11/06/2019 by Yosuke Tanigawa; add support for array_combined and array_imp_combined; code clean up.)
(Updated 11/27/2019 by Yosuke Tanigawa; add --keep option; code clean up.)

optional arguments:
  -h, --help            show this help message and exit
  --plink1              Flag to run GWAS with PLINK v1.9 instead of v2.0
  --run-array           Run GWAS on directly genotyped (array) data
  --run-exome           Run GWAS on exome data (Regeneron calls)
  --run-exome-gatk      Run GWAS on exome data (GATK calls)
  --run-imputed         Run GWAS on imputed data
  --run-cnv             Run GWAS on array-derived CNV genotypes
  --run-cnv-burden      Run CNV burden test (GWAS on 0/1 CNV overlaps gene,
                        from array-derived CNV genotypes)
  --run-hla             Run GWAS on imputed HLA allelotypes
  --run-array-imp-combined
                        Run GWAS on the array_imp_combined dataset
  --run-array-combined  Run GWAS on the array_combined dataset
  --pheno PHENO         Path to a phenotype file
  --out OUTDIR          Path to desired output *directory*. Summary stats will
                        be output according to phenotype name (derived from
                        passed file) and Rivas Lab specification for GBE.
                        Defaults to current working directory.
  --population POP      Flag to indicate which ethnic group to use for GWAS.
                        Must be one of all, white_british, non_british_white,
                        e_asian, s_asian, african
  --keep KEEP           A file that specifies the list of individuals for
                        GWAS. It can be overwritten by sex-specific
                        subcommands
  --remove-add RMADD    Flag to indicate file of list of additional
                        individuals to remove
  --keep-related        Flag to keep related individuals in GWAS. Default is
                        to remove them.
  --cores CORES         For underlying plink command/batch job submission:
                        Amount of cores to request. Default is 4 cores for
                        batch jobs.
  --memory MEM          For underlying plink command/batch job submission:
                        Amount of memory (in MB) to request. Default is 24000
                        for batch jobs.
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
  --sex-div SEX_DIV     Whether to run sex-divided GWAS, this removes sex as a
                        covariate from the analysis. Default is False. If
                        specifying this, please also use --keep-sex to specify
                        the sex to keep and --keep-sex-file to specify the
                        file location.
  --keep-sex KEEP_SEX   Which sex to keep, will be used in constructing the
                        output filenames, for use with --sex-div.
  --keep-sex-file KEEP_SEX_FILE
                        Location of the file specifying the IIDs to include
                        related to that sex, for use with --sex-div.
  --include-x INCLUDE_X
                        Whether to include the X chromosome, defaults to False
  --additional-plink-opts [PLINK_OPTS [PLINK_OPTS ...]]
                        Addtional options for plink
```
 Some examples are below:

- Input phenotype file: `--pheno INI50.phe`
- Output and logging directories: `--out ~/gwas/ --log-dir ~/gwas/logs/`
- Runtime options: `--run-array-combined`, `--run-imputed`, `--run-now`, `--run-exome`
- QC: `--population white-british --keep-related`

Our most common use cases of the script can be found in the [`run_array_hla_cnv.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/run_array_hla_cnv.sh) and [`run_exome.sh`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/run_exome.sh) wrapper scripts.

#### Notes: 

- Output and logging options must be directories; the script will automatically handle filename generations.
- The flag to keep related individuals currently only works with `--population all`, since the lab's population definitions are only on the set of unrelated individuals used for PCA.
- _Make sure to check the default parameters before running GWAS!_

### GWAS freeze

This part of the documentation is meant to take you through how to upload a "tarball" of summary statistics (which are capped at a certain user-specified p-value threshold) to the lab Google Drive.

1. Run the [check script](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/check_gwas/check_gwas.sh) or one of its wrappers (for [array](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/check_gwas/check_array_gwas.sh) or [exome](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/check_gwas/check_exome_gwas.sh)) and deposit the results to the [`check_gwas/output/`](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/check_gwas/output) subdirectory.
2. Inspect the results to make sure that there's nothing you don't expect (e.g., files with way too many NA lines, or `.phe` files that don't have summary statistics).
3. Run the [freeze script](https://github.com/rivas-lab/ukbb-tools/tree/master/04_gwas/gwas_freeze.sh), whose usage is described above, to generate a population-level "tarball" with proper ownership, and additionally generate a master PheWAS file for [`ukbb-query`](https://github.com/rivas-lab/ukbb-query).
4. The freeze should be available at `/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/freeze/YYYYMMDD/[pop]`.
5. Create an appropriate subdirectory within the lab Google Drive, like so: `rivas-lab/ukbb24983/cal/gwas/freeze/YYYYMMDD/[pop]/`.
5. Make sure the `rclone` module is configured as described [here](https://github.com/rivas-lab/wiki/wiki/Google-Drive-Storage#21-upload-with-rclone). Then, transfer and deposit the results (both the `.tsv.gz` and its corresponding tabix index file, `.tsv.gz.tbi`) to Google Drive using the instructions there.
6. Then, create a directory on Google Drive to store all of the sumstats files. An example naming convention would be `rivas-lab/ukbb24983/cal/gwas/freeze/20190913/e_asian/ukb24983_v2_hg19.e_asian.genotyped.glm.20190913`. Replace `20190913` with the current `YYYYMMDD`.
7. Finally, upload all of the sumstats files to the directory using the same check file you generated before. An example is below:

```
cat check_array_gwas.50195103_e_asian.out | awk '(NR>1){print $3}'| parallel --eta -j10 --header 1 rclone copy {} gdrive:rivas-lab/ukbb24983/cal/gwas/freeze/20190913/e_asian/ukb24983_v2_hg19.e_asian.genotyped.glm.20190913 :::: /dev/stdin
```

**DEPRECATED - Alternate `gdrive` module usage for upload:**

```
cat check_array_gwas.50195103_e_asian.out | awk '(NR>1){print $3}'| parallel --eta -j10 --header 1 gdrive upload -p 1M2tIqCavPZoNA3K4nX3rRB_sC3YB4vND {} :::: /dev/stdin
```

8. Finally, when relevant, transfer the files to GBE server(s) and load them to the SciDB database(s).