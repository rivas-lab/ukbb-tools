# UKB-wide genetic correlation (rg) with LDSC
## Yosuke Tanigawa, 2020/9/22

We use LD score regression (LDSC) to compute the genetic correlations among UK Biobank traits.

Here, we have [`1_ldsc_rg.sh`](1_ldsc_rg.sh), a wrapper that takes a pair of `GBE_ID`s and (optional argument of) population label (by default, we use the meta-analyzed summary statistics). Note that our LD score reference is computed for the European samples, which means the genetic correlation estimates for non-European samples may not be reliable.

## Data location

We used the meta-analyzed summary statistics and computed the LDSC rg.
The lateset intermediate results were written to 

`/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/rg.metal.20201015-115720.tsv` (7,618,633 files)

## Wrapper scripts

### Usage of the wrapper script

To run LDSC `--rg` for a pair of traits, please call the wrapper script with two GBE IDs.

```{bash}
bash 1_ldsc_rg.sh INI21001 HC221
```

By default, we use the meta-analyzed summary statistics. If you want to focus on particular population, please add an optional 3rd argument.

```{bash}
bash 1_ldsc_rg.sh INI21001 HC221 white_british
```

The results from this scripts will be written to: `/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/rg/<population>.<GBE_ID_1>.<GBE_ID_2>.log`

For example, the two example commands above results in the following two files:

- `metal.INI21001.HC221.log`
- `white_british.INI21001.HC221.log`

### Batch computation

To apply the script in scale, we first create the list of phenotypes and the non-redundant set of phenotype pairs.
To generate the list of phenotypes, we use [`2a_generate.GBE_ID.lst.sh`](2a_generate.GBE_ID.lst.sh). From this, we use `itertools.combinations(lst, 2)` in Python, where `lst` is the list of phenotypes generated from the previous script.

The resulting files were not commited to the GitHub repository because of the file size.

We then have a wrapper script that takes an index (of an array job) and optional `offset` argument to specify the starting index.

```{bash}
bash 2b_wrapper.sh <array job index> <offset>
```

You can use SLURM jobs to compute the catalog of LDSC rg.

### Tabulate the results

We have another wrapper to tabulate the results.

```{bash}
bash 3_tabulate.sh
```

This script shows the results into a table format. When you finish your computation, please run this script and save the results into a file.

The latest `rg` results file can be found at `/oak/stanford/groups/mrivas/ukbb24983/array-combined/ldsc/rg.metal.20201102-173634.tsv`.

