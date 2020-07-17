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

## Move the results


cat 2_gwas.jobs.status.20200713-104853.tsv |  awk -v FS='\t' -v OFS='\t' -v T="TRUE" 'NR==1 || ($5 == T && $6 == T){print $8}' | while read f ; do if [ -f $f ] ; then echo $(zcat $f | wc -l) $f ; fi ; done

```
find data -name "*.gz" | while read f ; do echo $(zcat $f | wc -l) $f ; done
```