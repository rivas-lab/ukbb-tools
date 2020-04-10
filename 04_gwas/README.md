# Genome-wide association studies with UK Biobank data

The "main" script in this directory, `gwas.py`, is designed to be a convenient wrapper for running genome-wide associations with PLINK using the lab's resources from the UK Biobank. This directory contains a plethora of options and double-wrapper scripts around `gwas.py`.

## Contents

1. `04_gwas_misc.sh`: A file containing common functions used across the GWAS scripts. Loaded by running `source 04_gwas_misc.sh`.
2. `flipfix_and_liftOver.py`: Fixes "flips" in coordinates and (optionally) lifts them over to hg38. *Known issue: Multiallelic sites are not supported.*
- Inputs: 
- Outputs: 
- Example usage:
3. `gwas.py`: The bread and butter of this directory. This script uses a Python command-line interface to create the recommended PLINK script and submits it to the cluster.
- Inputs: 
- Outputs: 
- Example usage:
4. `gwas_freeze.sh`: 
- Inputs:
- Outputs:
- Example usage:
5. `plink_updater.sh`: 
- Inputs:
- Outputs:
- Example usage: 
6. `run_array_hla_cnv.sh`: The default workhorse for running GWAS across phenotypes; as indicated by the name, it runs GWAS across the array, HLA, and CNV data.
- Inputs:
- Outputs:
- Example usage:
7. `run_exome.sh`: 
- Inputs:
- Outputs:
- Example usage:
8. `show_gwas.sh`: 
- Inputs:
- Outputs:
- Example usage:
9. `test-flipfix_and_liftOver.sh`: 
- Inputs:
- Outputs:
- Example usage:
10. `check_gwas`: 
11. `extras`: 
12. `flipfix`: 
13. `mvp`: 
14. `ref_alt`: 
15. `rerun_logs`: This directory contains the log files from the SLURM submissions of the `run_gwas.sh`-type files. None of the logs will show up on the GitHub due to a nifty `.gitignore`.

## Pipelines and Workflows

### GWAS options

A full list of options can be obtained by running `python gwas.py -h`, and some examples are below:

- Input phenotype file: `--pheno INI50.phe`
- Output and logging directories: `--out ~/gwas/ --log-dir ~/gwas/logs/`
- Runtime options: `--run-array`, `--run-imputed`, `--run-now`
- QC: `--population white-british --keep-related`
- Job submission options: `--sb-mem 16000 --sb-time 1-00:00:00 --sb-parti normal,mrivas`

### Notes: 

- Output and logging options must be directories — the script will automatically handle filename generations
- Either or both of `--run-array`, and `--run-imputed` can be specified, but `--run-now` (which runs the analysis immediately rather than submitting as a batch job) is only compatible with `--run-array`.
- The flag to keep related individuals currently only works with `--population all`, since the lab's population definitions are only on the set of unrelated individuals used for PCA.
- _Make sure to check the default parameters before running!_

### GWAS freeze

Please ask the core lab members for this freeze procedure. 
Here is a brief summary of the procedure.

##### 1) run the check script and deposit the results to `check_gwas/output/` subdir

##### 2) inspect the results

##### 3) run the freeze script: 

```
$ # you can check the usage:
$ bash gwas_freeze.sh 
$ # example
$ bash gwas_freeze.sh check_gwas/output/check_array_gwas.50195103_e_asian.out ukb24983_v2_hg19 genotyped 20190913 e_asian
```

This creates the tar ball with proper ownership and generate master PheWAS file for `ukbb-query`.
https://github.com/rivas-lab/ukbb-query

##### 4) check the results. The freeze should be availble at 

Check the number of lines, number of non-NA lines, etc.

`/oak/stanford/groups/mrivas/ukbb24983/cal/gwas/freeze/20190913/e_asian`

##### 5) Transfer and deposit the results to google drive

https://github.com/rivas-lab/wiki/wiki/Google-Drive-Storage

Please create the relevant directory.

Example: `rivas-lab/ukbb24983/cal/gwas/freeze/20190913/e_asian/`

Please deposit the master sumstats file and its tabix index file to the directory you've just created (`rivas-lab/ukbb24983/cal/gwas/freeze/20190913/e_asian/`).

```
cd /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/freeze/20190913/e_asian
ml load gdrive
gdrive upload -p 1uzX4misJ-2q52ECmSzNwq8tPbMlqoc3Q check_array_gwas.50195103_e_asian.out
gdrive upload -p 1uzX4misJ-2q52ECmSzNwq8tPbMlqoc3Q ukb24983_v2_hg19.e_asian.genotyped.glm.20190913.1e3.tsv.gz.tbi
gdrive upload -p 1uzX4misJ-2q52ECmSzNwq8tPbMlqoc3Q ukb24983_v2_hg19.e_asian.genotyped.glm.20190913.1e3.tsv.gz
```

Then, create a directory to store each full sumstats file.

Example: `rivas-lab/ukbb24983/cal/gwas/freeze/20190913/e_asian/ukb24983_v2_hg19.e_asian.genotyped.glm.20190913`

Then upload all of the sumstats files in the directory:

```
cat check_array_gwas.50195103_e_asian.out | awk '(NR>1){print $3}'| parallel --eta -j10 --header 1 gdrive upload -p 1M2tIqCavPZoNA3K4nX3rRB_sC3YB4vND {} :::: /dev/stdin
```

##### 6) Trasnfer and load the files to GBE (when relevant)

When relevant, please transfer the files to GBE server(s) and load them to the SciDB DB(s).