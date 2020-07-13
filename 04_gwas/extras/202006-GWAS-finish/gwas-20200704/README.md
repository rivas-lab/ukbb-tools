# GWAS

## GWAS finalize job submission to generate 236 files.

```{bash}
cat ../gwas-current-gz-wc.20200704-155715.combined.tsv | awk 'NR == 1 || $5 < 1059397' | egrep -v 'INI21049|INI21051|INI21052|INI21053|INI21054|INI21055|INI21056|INI21058|INI21059|INI21060|INI21061' | cut -f1-2  > 1_gwas.jobs.tsv
```

```{bash}
sbatch --mail-user=ytanigaw@stanford.edu --mail-type=END,FAIL --account=mrivas -p nih_s10 --cpus-per-task=6 --ntasks=1 --mem=30000 --time=2-0:00:00 --job-name=gwas --output=logs/gwas.%A_%a.out --error=logs/gwas.%A_%a.err --array=1-236 1_gwas.sh

Submitted batch job 15922450

scontrol update job 15922450 partition=nih_s10,batch
```

##

{
    echo "#GBE_ID population" | tr ' ' '\t'
    cat 1_gwas.jobs.tsv | awk 'NR>1' | while read GBE_ID pop ; do
    f=data/${pop}/ukb24983_v2_hg19.${GBE_ID}.array-combined.log
    if [ -s $f ] ; then echo ${GBE_ID} ${pop} | tr ' ' '\t' ; fi
    done
} > 1_gwas.jobs.finished.$(date +%Y%m%d-%H%M%S).tsv
