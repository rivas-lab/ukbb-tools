# variant annotation

This is Yosuke's attempts to quickly generate variant annotation.
We performed this analysis to get some annotation for non-autosomal variants.

## file location

`/oak/stanford/groups/mrivas/ukbb24983/array_combined/annotation/ukb24983_cal_hla_cnv.non-autosomes.vep.tsv`

## How to convert pvar to VCF file

- [`pvar.to.vcf.sh`](pvar.to.vcf.sh): this script takes pvar[.zst|.gz] file and convert it to a VCF file.

### usage

```{bash}
bash pvar.to.vcf.sh /oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv.pvar.zst | awk '$4 !="N"' | head -n20  > data/ukb24983_cal_hla_cnv.head.vcf
```

## How to use VEP

### prepare reference data (one-time set up)

We have our reference data in `/oak/stanford/groups/mrivas/public_data/vep_cache_20170410`. Please make a symlink from a home directory.

```{bash}
cd $HOME
ln -s /oak/stanford/groups/mrivas/public_data/vep_cache_20170410 .vep
```

- Note that the reference data is too old.

### running the VEP command

```{bash}
ml load perl
perl \
/oak/stanford/groups/mrivas/users/guhan/software/ensembl-vep/vep.pl \
-offline --force_overwrite \
-i data/ukb24983_cal_hla_cnv.head.vcf \
-o data/ukb24983_cal_hla_cnv.head.vep.tsv
```

## Variant annotation

### non-autosomal variants in the array combined dataset

```{bash}
 bash pvar.to.vcf.sh /oak/stanford/groups/mrivas/ukbb24983/array_combined/pgen/ukb24983_cal_hla_cnv.pvar.zst | awk '$4 !="N"' | sed -e "s/chrXY/chrX/g" | egrep '^#|^chrX|^chrY|^chrMT' > data/ukb24983_cal_hla_cnv.non-autosomes.vcf

 perl /oak/stanford/groups/mrivas/users/guhan/software/ensembl-vep/vep.pl -offline --force_overwrite -i data/ukb24983_cal_hla_cnv.non-autosomes.vcf -o data/ukb24983_cal_hla_cnv.non-autosomes.vep.tsv
```

#### ToDo

- This non-autosomal variant annoation does not have Loftee etc.

## reference

- [Docker VEP](https://hub.docker.com/r/ensemblorg/ensembl-vep)
  - Perhaps, we may want to use the latest VEP via Docker/Singularity?
