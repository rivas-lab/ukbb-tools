# coding 339

- 11 phenotypes with coding 339 was annotated in a wrong way.
- INI21049,INI21051,INI21052,INI21053,INI21054,INI21055,INI21056,INI21058,INI21059,INI21060,INI21061

- `/oak/stanford/groups/mrivas/private_data/ukbb/24983/phenotypedata/2005693/37855/` is the table file

We will generate a new phe file for those.

- [`1_extract_coding339.R`](1_extract_coding339.R): generate a dump from a tab file
- [`2_phe_definition.ipynb`](2_phe_definition.ipynb): generate the phe files under `data/`
  - `data/ukb.2005693.37855.coding339.phe`: the combined phe file
  - We have single phe files:

```
data/ukb.2005693.37855.coding339.INI21048.phe
data/ukb.2005693.37855.coding339.INI21049.phe
data/ukb.2005693.37855.coding339.INI21051.phe
data/ukb.2005693.37855.coding339.INI21052.phe
data/ukb.2005693.37855.coding339.INI21053.phe
data/ukb.2005693.37855.coding339.INI21054.phe
data/ukb.2005693.37855.coding339.INI21055.phe
data/ukb.2005693.37855.coding339.INI21056.phe
data/ukb.2005693.37855.coding339.INI21058.phe
data/ukb.2005693.37855.coding339.INI21059.phe
data/ukb.2005693.37855.coding339.INI21060.phe
data/ukb.2005693.37855.coding339.INI21061.phe
```

- [`3_gwas.sh`](3_gwas.sh): run GWAS

Test with e_asian

```{bash}
bash 3_gwas.sh e_asian
```

Confirmed that the line numbers look good (1080969).

```{bash}
sbatch -p mrivas --qos=high_p --nodes=1 --mem=32000 --cores=6 --time=1-0:00:00 --job-name=gwas339 --output=logs/gwas339.%A.out --error=logs/gwas339.%A.err 3_gwas.sh white_british

Submitted batch job 3572796
```

```{bash}
for pop in 'related' ; do bash 3_gwas.sh ${pop} ; done
```

```{bash}
for pop in 'others' 'non_british_white' 'african' 's_asian' ; do sbatch -p mrivas --qos=high_p --nodes=1 --mem=32000 --cores=6 --time=1-0:00:00 --job-name=gwas339.${pop} --output=logs/gwas339.${pop}.%A.out --error=logs/gwas339.${pop}.%A.err 3_gwas.sh ${pop} ; done
```

Submitted batch job 3575963
Submitted batch job 3575964
Submitted batch job 3575965
Submitted batch job 3575966
