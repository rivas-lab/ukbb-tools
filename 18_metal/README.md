# Meta-analysis with Metal

[We have Metal installed on Sherlock](https://github.com/rivas-lab/sherlock-modules/tree/master/metal).

This directory contains helper scripts to run Metal.

## ToDo

- Add post-processing scripts so that metal file has CHROM and POS column.
- Apply flipfix
- Apply `bgzip`


## `run_metal.sh`

This is the helper script to run Metal.

### Usage of `run_metal.sh`

Let's say you have a list of summary statistics. You will pass the list (file) as the first argument of the script.

The script generates the masterfile (used as the input for Metal) and call Metal. Metal creates files starting with the specified prefix.

In the following example, `hearing.lst` contains the list of your sumstats, `hearing.masterfile` is the output from the script, and `$(pwd)/out/METAANALYSIS_hearing_` is the prefix of the output files from Metal.

```{bash}
bash run_metal.sh hearing.lst hearing.masterfile $(pwd)/out/METAANALYSIS_hearing_
```

## `generate_master_file.sh`

If you don't want to run Metal but just want to generate the masterfile, you may use this script.

### Usage of `generate_master_file.sh`

```{bash}
bash run_metal.sh hearing.lst hearing.masterfile
```

## `add_BETA_from_OR.sh`

For plink summary statistics for binary outcomes, we have OR column instead of BETA column.

This helper script adds BETA column by computing `log(OR)`.

### Usage of `add_BETA_from_OR.sh`

```{bash}
bash add_BETA_from_OR.sh BBJ.Hearing_Loss.gz BBJ.Hearing_Loss.BETA.gz
```
