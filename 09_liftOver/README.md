
```
cat ukb24983_cal_cALL_v2.bim \
| sed -e 's/^23/X/g' | sed -e 's/^24/Y/g' | sed -e 's/^25/X/g' | sed -e 's/^26/M/g' \
| awk -v OFS='\t' -v FS='\t' '{print "chr" $1,$4,$4+1,$2,$3,$5,$6}' \
| liftOver /dev/stdin /cluster/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz hg38 hg38.unmapped

python generate_liftOver_tsv.py
bgzip ukb24983_cal_cALL_v2.liftOver.tsv
```

