# GRM

Using GCTA, we compute GRM

https://cnsgenomics.com/software/gcta/#Overview


## job submission

```
[ytanigaw@sh-102-07 ~/repos/rivas-lab/ukbb-tools/15_GRM]$ sbatch --array 2-1000 gcta_grm.sbatch.sh
Submitted batch job 55066375
```

### combine

```
# Merge all the parts together (Linux, Mac)
cat test.part_3_*.grm.id > test.grm.id
cat test.part_3_*.grm.bin > test.grm.bin
cat test.part_3_*.grm.N.bin > test.grm.N.bin
```


## debug code

```
[ytanigaw@sh-102-07 ~/repos/rivas-lab/ukbb-tools/15_GRM]$ gcta64 --make-grm-part 22 20 --bfile /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/
ukb24983_cal_chr22_v2 --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --maf 0.01  --thre
ad-num 6 --out test.20
*******************************************************************
* Genome-wide Complex Trait Analysis (GCTA)
* version 1.92.4 beta Linux
* (C) 2010-2019, The University of Queensland
* Please report bugs to Jian Yang <jian.yang@uq.edu.au>
*******************************************************************
Analysis started at 21:53:35 UTC on Wed Nov 20 2019.
Hostname: sh-102-07.int

Options:

--make-grm-part 22 20
--bfile /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_chr22_v2
--keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe
--maf 0.01
--thread-num 6
--out test.20

```

