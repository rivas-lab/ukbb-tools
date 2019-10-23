# LD map

We compute some files representing linkage disequilibrium (LD) structure. 
Specifically, we generate the following for each population for each dataset type.

- `bool.prune.in`
- `bool.orune.out`
- `ld_map.ld.gz`

The first two files are the results of LD pruning. 

The last one is the output from `--r2`.

## some scripts in this directory

### `imp_v3_ldmap.sbatch.sh` 
This is a job script to compute LD map for imputation dataset

### `imp_v3_ldmap.check.sh`
This is a check script to see if we have the results

