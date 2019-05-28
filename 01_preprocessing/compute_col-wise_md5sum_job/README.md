
##  How the job file is generated

```
find /oak/stanford/groups/mrivas/ukbb/24983/phenotypedata/ -name "*.tab.columns" | sort | grep -v old | while read f ; do wc -l $f ; done | while read n_cols f ; do (( n_jobs = n_cols / 50 + 1 )) ; for i in $(seq 1 $n_jobs) ; do echo $i ${f%.columns} ; done  ; done > jobs.tsv
```

