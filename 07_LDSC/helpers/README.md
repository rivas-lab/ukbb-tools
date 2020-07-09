# LDSC helper scripts

In 2020/4, we rewrote wrapper scripts for LDSC. Previous scripts can be found in [`old_helpers`](old_helpers) directory.

If you're interested in using our helper script, please first read [LDSC documentation](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) so that you have some understanding on what scripts needs to be called.

We now use Docker/Singularity image of LDSC (rather than assuming everybody has their LDSC environment). The following three scirpts are basically wrappers for Dockerlized LDSC.

- [`ldsc_munge.sh`](ldsc_munge.sh)
  - This scripts internally calls [`make_ldsc_input_file_v2.R`](make_ldsc_input_file_v2.R) to filter the input summary statistics.
  - We use `--merge` for the imputed data, but not for the array data. Unless you are using the imputed dataset for LDSC, you'll likely to [lose a lot of variants](https://gist.github.com/yk-tanigawa/f7895b1b2bd94facee4774e3e39f3c14).
- [`ldsc_h2.sh`](ldsc_h2.sh)
- [`ldsc_rg.sh`](ldsc_rg.sh)

To see the usage of each script, please execute each script without arguments.

Please check [our documentation on the module file](https://github.com/rivas-lab/sherlock-modules/tree/master/ldsc) to understand the environmental variables used in the scripts.

To tabulate the results, we have

- [`ldsc_rg_view.sh`](ldsc_rg_view.sh)
