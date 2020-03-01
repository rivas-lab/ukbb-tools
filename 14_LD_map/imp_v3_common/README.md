# LD map for common variants in imputed dataset (v3)

## results

- `/scratch/groups/mrivas/ukbb24983/imp/ldmap_common`
  - This is symlinked from `$OAK`: `/oak/stanford/groups/mrivas/ukbb24983/imp/ldmap_common`

## scripts

### `imp_v3_ldmap.sbatch.sh`

This is a job script to compute LD map for imputation dataset. There are two variants of this script: `ldmap.lowmem.sbatch.sh` (low mem) and `ldmap.chr.sbatch.sh` (one for X and XY chr).

### `imp_v3_ldmap.check.sh`

This is a check script to see if we have the results

### `post_processing.sh`

This script was used to apply `bgzip` and `tabix`.
