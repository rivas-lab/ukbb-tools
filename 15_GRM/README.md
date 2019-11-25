# GRM

Using GCTA, we compute GRM

https://cnsgenomics.com/software/gcta/#Overview


## job submission

```
[ytanigaw@sh-102-07 ~/repos/rivas-lab/ukbb-tools/15_GRM]$ sbatch --array 2-1000 gcta_grm.sbatch.sh
Submitted batch job 55066375
```

### failed jobs

It seemed like 3/1,000 jobs failed (998,999,1000) after two trials.
It is likely that there is a bug in GCTA which causes seg fault.


## combine

```
gcta_grm.combine.sh
```



