# PheWAS

This directory contains scripts for running PheWAS and generating the files necessary to do so.

## Dependencies:

_Note_: this section only needs to be updated when new phenotype data comes in, but feel free to read it for more info about where the files used by `phewas.py` come from.

In order to run PheWAS, we need a master table containing all the phenotype data we've processed (see the [phenotyping](https://github.com/rivas-lab/ukbb-tools/blob/master/02_phenotyping/) folder for more info). To do this, we leverage the phenotype info table, which is also in the phenotyping folder. The maker script, `combine_phe.py` will walk through this file and selecting phenotypes with n/N > 100 for inclusion. It then loads each of those files into memory (!!) via a dictionary keyed on sample ID, then writes all the data out to file. 

If you're rerunning this, you'll probably want a lot of memory -- 64GB is sufficient in my experience (Matt).

### Master phenotype info table schema

The list of columns in the master phenotype info file is as follows:
```
'GBE_phe_code',
'GBE_phe_name',
'UKB_field ID',
'UKB_table ID',
'UKB_basket ID',
'UKB_application ID',
'N_total',
'N_White_British',
'N_African',
'N_East_Asian',
'N_South_Asian',
'tsv_file',
'version',
'phe_file'
```

## Analysis:

PheWAS can be run using `phewas.py`. All options can be viewed by running `phewas.py -h`. The script can be run in one of the following modes: gene, variant, region.

1. For variant mode simply provide a list of variant IDs after the `--variants` flag. These can be listed after the flag, separated by spaces, or you can pass a file with the variants you'd like to run. The script will tell you if any of the requested variants were not found in UK Biobank, but you may want to double-check that you've got the right IDs by looking them up in the `.bim` files first.

2. For region mode, specify a chromosomal window after the `--region` flag, formatted like so: `CHROM:BP1-BP2`. A gentle reminder that UK Biobank uses hg19 genomic coordinates.

3. For gene mode, specify the gene name after the `--gene` flag. The script will look up the corresponding region according to HGNC nomenclature and proceed accordingly.

_Important_: There is no limit as to the number of variants that can be specified by a PheWAS. As such, you'll need to be careful to allocate sufficient time and memory for the analysis. The script will warn you if your input parameters include over 100 variants, which is about the most I (Matt) would recommend. For reference, a PheWAS for _GDF15_ (45 variants, ~2500 phenotypes) took about 4 hours to run, using 16GB memory.