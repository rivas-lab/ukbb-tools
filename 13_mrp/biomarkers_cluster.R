library("tidyverse")
library("comprehenr")

#setwd("~/Downloads/")

corrs <- read.delim("biomarkers_array_rg.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
sumstats <- read.delim("sumstat_paths.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
phens <- sumstats$GBE_ID
R_phen <- matrix(0L, nrow = length(phens), ncol = length(phens))
print(phens)
for (row in 1:nrow(corrs)) {
  i = match(c(corrs[row, "p1"]), phens)
  j = match(c(corrs[row, "p2"]), phens)
  R_phen[i, j] = abs(corrs[row, "rg"])
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
plot(dist_hclust)
rect_hclust <- rect.hclust(dist_hclust, border="red")

count = 0
list1 = c()
list2 = c()
for (i in 1:length(rect_hclust)) {
  for (j in 1:length(rect_hclust[[i]])) {
    count = count + 1
    list1[count] <- rect_hclust[[i]][j]
    list2[count] <- i
  }
}
pearson_nohla_clusters <- data.frame("PHEN" = to_vec(for(i in 1:length(list1)) corrs_phens[list1[i]]), cluster = list2) %>% inner_join(phenos_to_include)

clusters <- pearson_nohla_clusters %>% inner_join(phenos_to_include)

write.table(clusters, file='clusters.tsv', quote=FALSE, sep='\t', row.names=FALSE)

cluster_sizes <-clusters %>% group_by(cluster) %>% count()
cluster_sizes <- cluster_sizes %>% filter(n >= 3 & n <= 50)
write.table(cluster_sizes, file='clusters_to_run.tsv', quote=FALSE, sep='\t', row.names=FALSE)

corrs_nohla_sig <- corrs %>% filter(PEARSON_P_NOHLA <= 5*10^-8) %>% select(PHEN_1, shortname_1, PHEN_2, shortname_2, PEARSON_NOHLA, PEARSON_P_NOHLA, BETALEN_NOHLA)
write.table(corrs_nohla_sig, file='corrs_pearson_nohla_sig.tsv', quote=FALSE, sep='\t', row.names=FALSE)
"""
