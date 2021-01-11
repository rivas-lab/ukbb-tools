# Multiple Rare-variants and Phenotypes Mixed Model

## Contents

1. [`MVMetaLearnMRP2016OctoberMRPMM.pdf`]():
2. [`gen_var_list.py`]():
- Inputs:
- Outputs:
- Example usage:
3. [`mrpmm.py`]():
- Inputs:
- Outputs:
- Usage: 
- Example usage:
4. [`old_mrpmm.py`]():
- Inputs:
- Outputs:
- Example usage:
5. [`plotting.py`]():
- Inputs:
- Outputs:
- Example usage:

## Pipelines and Workflows

In order to run MRPMM, you need the following:

- A `--variants` argument, which is the path to a file containing the list of 
variants to include, one line per variant. This should have a header of "V".
- A `--phenotypes` argument, which is the path to a tab-separated file containing 
a list of summary statistic file paths, phenotype names, and whether or not to 
use the file in R_phen generation. An example is below:

```       
path        pheno        R_phen
/path/to/file1   pheno1     TRUE
/path/to/file2   pheno2     FALSE
```

- A `-metadata_path` argument containing a file that annotates variants with 
gene symbols, consequences, and HGVSp annotations. An example is below:

```       
V       gene_symbol     most_severe_consequence HGVSp  
1:69081:G:C     OR4F5   5_prime_UTR_variant     ""
```

- An optional. `--C` argument, which accepts one or more positive integers 
space-separated, indicating the number(s) of clusters with which to run MRPMM.
Default 3.

- An optional, `--se_thresh` argument, which caps the standard error for analyzed 
variants at a specified float. Default 100.

- An optional, `--maf_thresh` argument, which caps the minor allele frequency for 
analyzed variants at a specified float. Default 0.01.
