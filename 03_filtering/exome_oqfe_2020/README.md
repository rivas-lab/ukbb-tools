# The exome 200k dataset

Yosuke Tanigawa, 2020/10/26

## file location

`/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020`

- `download`: (symlink to `/scratch`)
- `ukb24983_exomeOQFE.{pgen,pvar.zst,psam}`

## methods summary

- [`1_download.sh`](1_download.sh): wrapper script for `gfetch`, submitted with [`1_download.sbatch.sh`](1_download.sbatch.sh)
- [`2_merge.sh`](2_merge.sh): we combine the per-chromosome genotype files into one file.

### preparing the bim file

```
wget -nd  biobank.ctsu.ox.ac.uk/crystal/crystal/auxdata/UKBexomeOQFEbim.zip
```

### `gfetch` software

```
wget -nd  biobank.ctsu.ox.ac.uk/showcase/util/gfetch
chmor 770 gfetch
```

