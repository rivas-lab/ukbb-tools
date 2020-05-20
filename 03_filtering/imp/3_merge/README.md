# Merge array + imputation data

We generate the `array_imp_combined` dataset in the following procedure.

## Data availability

`/oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen`

### file name prefix error fix

We fixed some errors in file name prefix:

```
find /oak/stanford/groups/mrivas/ukbb24983/array_imp_combined -name "*ukb24983_ukb24983*" -type f | tee files.txt
paste files.txt <(cat files.txt | sed -e "s/ukb24983_ukb24983/ukb24983_hg19/g") | while read f g ; do mv $f $g ; done
```

## `3-1_merge.sh`

### version 1

We merge the following datasets:

- 5,182,706 variants on imputed dataset (v3)
  - MAF >= 0.01
  - Imputation quality >= 0.7
  - Biallelic only
  - The variant position is not present on genotyping array dataset
  - Missingness 
  - HWE 1e-7
- Array genotypes (v2)
- HLA
- CNV

The memory footprint of PLINK 1.9 `--merge` is proportional to the length of variant IDs.

```
From: Chris Chang <chrchang523@gmail.com>
Date: Wed, Sep 18, 2019 at 10:10 AM
Subject: Re: Out of memory error in PLINK --merge
To: Yosuke Tanigawa <ytanigaw@stanford.edu>


The memory requirement for PLINK 1.9's variant ID array scales as <length of LONGEST variant ID> * <number of variants>.
If you don't have another way to perform the merge, a workaround is to temporarily assign shorter unique IDs, and then set the IDs back to what they were after you're done using PLINK 1.9 on the dataset.

On Wed, Sep 18, 2019 at 9:30 AM Yosuke Tanigawa <ytanigaw@stanford.edu> wrote:
Hi Chris,
I'm wondering if you have some way to estimate the memory requirement forplink1.9 --merge operation.
I'm trying to generate a PLINK bed file for 488k individuals and ~10M variants.I gave 300GB of mem and 10 threads but got the following error:
Do you happen to know if there is a way to estimate the mem requirements?Related question: does this affected by having long variant IDs? If so, what's the memory footprint if all variant IDs are < 80 chars?
==================================================
Warning: Unusually long variant ID(s) present.  PLINK 1.9 does not scale well
to length-80+ variant IDs; consider using a different naming scheme for long
indels and the like.
/var/spool/slurmd/job50574351/slurm_script: line 33: 430069 Killed                  plink ${plink_opts} --merge-list ${tmp_merge_list} --make-bed --keep-allele-order --out ${dir}/ukb24983_imp_cALL_v3_maf1
slurmstepd: error: Detected 1 oom-kill event(s) in step 50574351.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.
==================================================
Thank you very much.
Sincerely yours,Yosuke
```

We submitted this job as

```
sbatch 3-1_merge.sh merge.lst.v1.imp.w.HWE.filter.tsv /oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined/pgen/ukb24983_ukb24983_cal_hla_cnv_imp
```

### version 2

We revised the variant QC criteria and removed the HWE filter.

```
cat merge.lst.v1.imp.w.HWE.filter.tsv | sed -e "s/_QC/_QC_noHWE/g" > merge.lst.v2.imp.wo.HWE.filter.tsv
```

```
sbatch 3-1_merge.sh merge.lst.v2.imp.wo.HWE.filter.tsv /oak/stanford/groups/mrivas/ukbb/24983/array_imp_combined/pgen_v2/ukb24983_ukb24983_cal_hla_cnv_imp
```


## 3-2 `3-2_recover_names.ipynb`

We recover the variant IDs using the codes in the notebook.

## 3-3 `3-3_bed_to_pgen.sh`

```
sbatch 3-3_bed_to_pgen.sh
```
