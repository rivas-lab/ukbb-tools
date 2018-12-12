## Genome-wide associations with UK Biobank data

The script in this directory, `gwas.py`, is designed to be a convenient PLINK wrapper for running genome-wide associations with lab resources from the UK Biobank. 

A full list of options can be obtained by running `python gwas.py -h`, and a summary is below:

- File input: phenotype
- File output: logging and out directories
- Data options: Genotyped/Imputed data
- QC: population specification, relatedness filter (currently non-functional in conjunction with populations)
- Run now 
- batch job submission parameters (memory, time, partitions)
