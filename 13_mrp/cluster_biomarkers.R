library("tidyverse")
library("comprehenr")

setwd("~/Downloads/")

corrs <- read.delim("corrs_names.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE) %>% 
  inner_join(read.delim("betalens_names.tsv", sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE))
corrs <- corrs %>% filter(PHEN_1 != PHEN_2)
final_phenos <- read.delim("final_phe_codes.txt", sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
phenos_to_include <- corrs %>% filter(PHEN_1 %in% final_phenos$V1) %>% select(PHEN_1, shortname_1) %>% unique() %>% rename(PHEN = PHEN_1, shortname = shortname_1)
corrs <- corrs %>% filter(PHEN_1 %in% phenos_to_include$PHEN) %>% filter(PHEN_2 %in% phenos_to_include$PHEN) %>% filter(BETALEN_NOHLA >= 10)

pearson_nohla_corrs_sig <- corrs %>% filter(PEARSON_P_NOHLA <= 5*10^-8)
pearson_nohla_corrs_sig_phens <- unique(c(pearson_nohla_corrs_sig$PHEN_1, pearson_nohla_corrs_sig$PHEN_2))
pearson_nohla_corrs_sig <- pearson_nohla_corrs_sig %>% filter(PHEN_1 %in% pearson_nohla_corrs_sig_phens) %>% filter(PHEN_2 %in% pearson_nohla_corrs_sig_phens)
pearson_nohla_R_phen <- matrix(0L, nrow = length(pearson_nohla_corrs_sig_phens), ncol = length(pearson_nohla_corrs_sig_phens))

for (row in 1:nrow(pearson_nohla_corrs_sig)) {
  i = match(c(pearson_nohla_corrs_sig[row, "PHEN_1"]), pearson_nohla_corrs_sig_phens)
  j = match(c(pearson_nohla_corrs_sig[row, "PHEN_2"]), pearson_nohla_corrs_sig_phens)
  pearson_nohla_R_phen[i, j] = pearson_nohla_corrs_sig[row, "PEARSON"]
}

for (row in 1:length(pearson_nohla_corrs_sig_phens)) {
  for (col in 1:length(pearson_nohla_corrs_sig_phens)) {
    if (row > col) {
      pearson_nohla_R_phen[row, col] = pearson_nohla_R_phen[col, row]
    }
  }
}

pearson_nohla_distances <- dist(pearson_nohla_R_phen)
dist_hclust <- hclust(pearson_nohla_distances, method = "ward.D2")
plot(dist_hclust)
rect_hclust <- rect.hclust(dist_hclust, k = 100, border="red")
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
pearson_nohla_clusters <- data.frame("PHEN" = to_vec(for(i in 1:length(list1)) pearson_nohla_corrs_sig_phens[list1[i]]), cluster = list2) %>% inner_join(phenos_to_include)

clusters <- pearson_nohla_clusters %>% inner_join(phenos_to_include)

write.table(clusters, file='clusters.tsv', quote=FALSE, sep='\t', row.names=FALSE)

cluster_sizes <-clusters %>% group_by(cluster) %>% count()
cluster_sizes <- cluster_sizes %>% filter(n >= 3 & n <= 50)
write.table(cluster_sizes, file='clusters_to_run.tsv', quote=FALSE, sep='\t', row.names=FALSE)

corrs_nohla_sig <- corrs %>% filter(PEARSON_P_NOHLA <= 5*10^-8) %>% select(PHEN_1, shortname_1, PHEN_2, shortname_2, PEARSON_NOHLA, PEARSON_P_NOHLA, BETALEN_NOHLA)
write.table(corrs_nohla_sig, file='corrs_pearson_nohla_sig.tsv', quote=FALSE, sep='\t', row.names=FALSE)