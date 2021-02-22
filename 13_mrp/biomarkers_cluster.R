library("tidyverse")
library("comprehenr")
library("dynamicTreeCut")
library("reshape2")

setwd("~/Dropbox/ukbb-tools/13_mrp")

corrs <- read.delim("biomarkers_array_rg.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
sumstats <- read.delim("sumstat_paths.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
phens <- sumstats$GBE_ID
R_phen <- matrix(0L, nrow = length(phens), ncol = length(phens))
print(phens)
for (row in 1:nrow(corrs)) {
  i = match(c(corrs[row, "p1"]), phens)
  j = match(c(corrs[row, "p2"]), phens)
  R_phen[i, j] = corrs[row, "rg"]
}

for (row in 1:length(phens)) {
  for (col in 1:length(phens)) {
    if (row > col) {
      R_phen[row, col] = R_phen[col, row]
    } else if (row == col) {
      R_phen[row, col] = 1
    }
  }
}

distances <- dist(R_phen)
dist_hclust <- hclust(distances, method = "ward.D2")
plot(dist_hclust, labels=sumstats$TRAIT)
assignments <- cutreeDynamic(dist_hclust, method="tree", deepSplit = 4, minClusterSize = 5)
df <- data.frame(GBE_ID=phens, phenotype=sumstats$X.Phenotype, cluster=assignments)
df <- df %>% arrange(GBE_ID) %>% arrange(cluster)
write.table(df, file='biomarkers_clusters.tsv', quote=FALSE, sep='\t', row.names=FALSE)
