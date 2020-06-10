# Flip check of the resulting file.

To ensure the REF allele in the resulting merged file matches with the REF in the FASTA file, 
we applied the flip-check script.

## flipcheck for v2

[`flipcheck.v2.sh`](flipcheck.v2.sh) is the flipcheck script.

[`flipcheck.v2.out`](flipcheck.v2.out) is the output (where we see a mismatch between the REF allele in pvar file and the REF from the fasta file)

We did not check the ref/alt flip for HLA and CNVs.

## flipcheck for v1

```
ml load bedtools
bash flipcheck.sh > flipcheck.out 2> flipcheck.err
```

Note that we excluded the CNV alleles from the flipcheck (the second `awk` command in the script removing ``$4 == N` and `$4 == P`)


We found the resulting file is almost perfect with 5 exceptions (for autosomes). It seemed like we have some flips for chrM but it should be fine....

```
[ytanigaw@sh-109-53 /oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-tools/03_filtering/imp/5_ref_alt]$ cat flipcheck.out
#CHROM  POS     ID      REF     ALT     FASTA_REF
4       7435305 Affx-24919403   T       C       A
7       150896003       Affx-29977087   T       G       C
12      53011955        Affx-8165091    T       C       G
21      46354975        Affx-18140622   A       G       C
22      39218764        Affx-19589717   T       C       G
```
