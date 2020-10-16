# Meta-analysis with Metal for the Exome dataset
## Yosuke Tanigawa, 2020/10/14

We have metal wrapper configured for the exome dataset.

```{bash}
bash 1_metal.sh INI50
bash 1_metal.sh FH1065
```

The results are currently written to `/oak/stanford/groups/mrivas/ukbb24983/exome/metal/20201014`

Note: we assume the input files are available under `/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/<population>` where, population is one of the followings: `white_british`, `non_british_white`, `african`, `s_asian`, `e_asian`, `related`, and `others`.

The actual set of sumstats files used in the M-A is written in the `metal.info.txt` file.

