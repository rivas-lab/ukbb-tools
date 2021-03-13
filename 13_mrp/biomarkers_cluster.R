library("tidyverse")
library("dynamicTreeCut")
library("corrplot")

setwd("~/Dropbox/ukbb-tools/13_mrp")

corrs <- read.delim("biomarkers_array_rg.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
sumstats <- read.delim("sumstat_paths.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
corrs <- corrs %>%
  separate(p1, c("prefix","phen_1","array-combined","sumstats","gz"), sep="\\.") %>%
  select(phen_1, p2, rg, p) %>% separate(p2, c("prefix","phen_2","array-combined","sumstats","gz"), sep="\\.") %>%
  select(phen_1, phen_2, rg, p) %>% rename(p1 = phen_1, p2 = phen_2)
phens <- sumstats$GBE_ID
R_phen <- matrix(0L, nrow = length(phens), ncol = length(phens))

for (row in 1:nrow(corrs)) {
  i = match(c(corrs[row, "p1"]), phens)
  j = match(c(corrs[row, "p2"]), phens)
  if (corrs[row, "p"] <= 0.01) {
    R_phen[i, j] = corrs[row, "rg"]
  } else {
    R_phen[i, j] = 0
  }
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

colnames(R_phen) <- sumstats$TRAIT
row.names(R_phen) <- sumstats$X.Phenotype
corrplot(R_phen, method="square", order="hclust", tl.col = "black")

distances <- dist(R_phen)
dist_hclust <- hclust(distances, method = "ward.D2")
plot(dist_hclust, labels=sumstats$TRAIT)
assignments <- cutreeDynamic(dist_hclust, method="tree", deepSplit = 4, minClusterSize = 5)
df <- data.frame(GBE_ID=phens, phenotype=sumstats$X.Phenotype, cluster=assignments, abbrev=sumstats$TRAIT)

# Corrplots of rg-based clusters
clusters <- unique(df$cluster)

for (i in 1:length(clusters)) {
  traits <- (df %>% filter(cluster == toString(clusters[i])))$abbrev
  traits <- sumstats$TRAIT[sumstats$TRAIT %in% traits]
  trait_names <- (df %>% filter(cluster == toString(clusters[i])))$phenotype
  trait_names <- sumstats$X.Phenotype[sumstats$X.Phenotype %in% trait_names]
  corrplot(R_phen[trait_names, traits], method="square", order="hclust", tl.col = "black")
}

df <- df %>% arrange(GBE_ID) %>% arrange(cluster)
write.table(df, file='biomarkers_clusters.tsv', quote=FALSE, sep='\t', row.names=FALSE)
