## Genome-wide associations with UK Biobank data

The script in this directory, `gwas.py`, is designed to be a convenient wrapper for running genome-wide associations with PLINK, using the lab's resources from the UK Biobank. 

### Analysis ptions:
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


