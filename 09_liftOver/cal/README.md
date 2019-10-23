# liftOver hg19 array genotype call data to hg38


## Step 1. Apply UCSC liftOver

Input file is `$OAK/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2_hg19.bim`

We run the following commands:

```
{ echo "#CHR_hg38 POS_hg38 ID REF_hg19 ALT_hg19 CHR_hg19 POS_hg19" | tr " " "\t" ;
cat ukb24983_cal_cALL_v2_hg19.bim \
| awk -v OFS='\t' -v chr="chr" -v sep=":" \
'{print chr $1, $4-1, $4-1+length($5), $2 sep chr $1 sep $4 sep $5 sep $6}' \
| liftOver /dev/stdin \
/cluster/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz \
/dev/stdout \
ukb24983_cal_cALL_v2_hg19.liftover.hg38.unmapped.txt \
| tr ":" "\t" | awk -v FS='\t' -v OFS='\t' '{print $1, $2+1, $4, $8, $7, $5, $6}'
} | gzip -9 > ukb24983_cal_cALL_v2_hg19.liftover.hg38.flip.unchecked.tsv.gz

gzip -9 ukb24983_cal_cALL_v2_hg19.liftover.hg38.unmapped.txt
```

## Step 2. Check REF allele from FASTA file

```
zcat ukb24983_cal_cALL_v2_hg19.liftover.hg38.flip.unchecked.tsv.gz \
| sed -e "s/chr//g" \
| bash ../04_gwas/flipfix/flipcheck.sh /dev/stdin \
/oak/stanford/groups/mrivas/public_data/genomes/hg38/hg38.fa \
| sed -e "s/FASTA_REF/FASTA_hg38_REF/g" \
| gzip -9 > ukb24983_cal_cALL_v2_hg19.liftover.hg38.flip.checked.tsv.gz
```

## Step 3. Prepare hg38 PLINK files

See `ukb24983_cal_cALL_v2_hg38.ipynb`

## Step 4. Deposit the resulting files

We remove the intermediate files

```
rm ukb24983_cal_cALL_v2_hg19.liftover.hg38.flip.unchecked.tsv.gz
rm ukb24983_cal_cALL_v2_hg38_tmp*
rm ukb24983_cal_cALL_v2_hg38.log # the contents of the log file is in ipynb in step 3.
echo https://github.com/rivas-lab/ukbb-tools/blob/master/09_liftOver/ > ukb24983_cal_cALL_v2_hg38.log
```

We deposit the following files to `$OAK/ukbb24983/cal/pgen`
- `ukb24983_cal_cALL_v2_hg19.liftover.hg38.flip.checked.tsv.gz`
- `ukb24983_cal_cALL_v2_hg19.liftover.hg38.unmapped.txt.gz`
- `ukb24983_cal_cALL_v2_hg38.pgen`
- `ukb24983_cal_cALL_v2_hg38.pvar`
- `ukb24983_cal_cALL_v2_hg38.psam`
- `ukb24983_cal_cALL_v2_hg38.log`
