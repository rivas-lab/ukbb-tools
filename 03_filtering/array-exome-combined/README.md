## version history

### 20210111

- The current version.
- We properly handle the hg19/hg38 issue.
- We incoporate the LD map / LD indep results computed [elsewhere](/14_LD_map/array-exome-combined).

### 20201217

- This was our original attempt
- It turned out that there was an critical error in "plink1.9 --bmerge" operation. We passed the array-combined and the exome 200k genotype data matrices, but forgot to supply hg19 version of bim file for exome 200k. Because of this, the resulting (sorted) genotype file had wrong coordinates for most variants. Because there were so many analyses depends on this wrong genotype file (including the LD map computation), we decided to re-run the analysis.

