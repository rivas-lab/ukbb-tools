# Merge array + imputation data

We generate the `array_imp_combined` dataset in the following procedure.

## step 0. Start with biallelic SNVs w/ MAF>1% on the imputation dataset

## step 1. `3-1_compare_cal_and_imp.ipynb`

We compared the imputation data and array genotyped data. 
We identified `673,909` variants have the same coordinates based on `inner_join` on (CHR, POS) pair.
(Note: `imp` has 10,231,518 variants and `cal` (array) has 805,426 variants)

We saved the results of this join to `/oak/stanford/groups/mrivas/ukbb/24983/imp/pgen/maf1/ukb24983_cal_v2_hg19_imp_v3_maf1.join.tsv`

## `3-2_make_imp_maf1_no_cal.sh`

Using the results from the previous step, we generate PLINK 1.9 & PLINK 2 files for the imputation variants that are not part of the array (`cal`).

## `3-3_merge.sh`

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

The longest variant IDs are:
- cal 13
- hla 9
- cnv 24
- imp 677

We therefore assigned temporal variant IDs for the merge.

`3-3_merge.sh` calls `assign_short_names_for_merge.R` and perform merge on dataset with short IDs.


We submitted this job as

```
3-3_merge.sh merge.lst.tsv
```

## 3-4 `3-4_recover_names.ipynb`

```
[ytanigaw@sh-ln06 login /oak/stanford/groups/mrivas/ukbb24983/array_imp_combined/pgen]$ mv ukb24983_ukb24983_cal_hla_cnv_imp.bim ukb24983_ukb24983_cal_hla_cnv_imp.shortnames.bim
```

We recover the variant IDs using the codes in the notebook.

## 3-5 `3-5_variant_qc.sh`

We compute the basic statistics (allele frequency and HWE p-value) for this combined set.
For array variant, we should use allele frequncy that is aware of the batch info.
