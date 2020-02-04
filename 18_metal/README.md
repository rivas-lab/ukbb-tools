# Meta-analysis with Metal

[We have Metal installed on Sherlock as a module](https://github.com/rivas-lab/sherlock-modules/tree/master/metal).

This directory contains helper scripts to run Metal.

## `run_metal.sh`

This is the helper script to run Metal. Please type `bash run_metal.sh` to see the usage of the wrapper script.

```{bash}
run_metal.sh: 2 positional arguments are required
run_metal.sh (version 0.0.1)
Run run_metal

Usage: run_metal.sh [options] in_file_list outfile_prefix
  in_file_list      A file that has a list of input files for METAL
  metal_out_prefix  The prefix of output files from METAL

Options:
  --flipcheck_sh     The location of flip check script
  --nCores     (-t)  Number of CPU cores

Default configurations:
  flipcheck_sh=/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/09_liftOver/flipcheck.sh
  nCores=4
```

This script depends on flipcheck script (which should be pointed to the one under `09_liftOver` by default) to check and fix the allele flips in the output files from Metal.

### Usage of `run_metal.sh`

Let's say you have a list of summary statistics. You will pass the list (file) as the first argument of the script.

The script generates the masterfile (used as the input for Metal) and call Metal. Metal creates files starting with the specified prefix.

In the following example, `LpA.lst` contains the list of your sumstats and `meta.LpA` is the prefix of the output files from Metal.

```{bash}
bash run_metal.sh LpA.lst meta.LpA
```

The first argument should contain the list of summary statistics (they will be passed to Metal in the specified order). For example, the contents of `LpA.lst` is as follows:

```
/oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/w_british/ukb24983_v2_hg19.Lipoprotein_A.genotyped.glm.linear.gz
/oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/non_british_white/ukb24983_v2_hg19.Lipoprotein_A.genotyped.glm.linear.gz
/oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/african/ukb24983_v2_hg19.Lipoprotein_A.genotyped.glm.linear.gz
/oak/stanford/groups/mrivas/projects/biomarkers_rivas/main/s_asian/ukb24983_v2_hg19.Lipoprotein_A.genotyped.glm.linear.gz
```

## `add_BETA_from_OR.sh`

For plink summary statistics for binary outcomes, we have OR column instead of BETA column.

This helper script adds BETA column by computing `log(OR)`.

### Usage of `add_BETA_from_OR.sh`

```{bash}
bash add_BETA_from_OR.sh BBJ.Hearing_Loss.gz BBJ.Hearing_Loss.BETA.gz
```

## `generate_master_file.sh` (for more customized analysis with Metal)

While `run_metal.sh` provides an automated analysis, one can manually edit the input file (`mastefile`) for Metal.
This script (`generate_master_file.sh`) generates a master file for the list of input summary statistic files.

### Usage of `generate_master_file.sh`

```{bash}
bash generate_master_file.sh LpA.lst meta.LpA
```
