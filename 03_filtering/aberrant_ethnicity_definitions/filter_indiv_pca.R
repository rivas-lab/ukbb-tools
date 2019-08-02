sqc.header <- read.table("/oak/stanford/groups/mrivas/ukbb24983/sqc/download/ukb_sqc_v2.fields.txt", header=F, stringsAsFactors=F)$V1
fam <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2.fam", header=F)
sqc <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/sqc/download/ukb_sqc_v2.txt", header=F)
colnames(sqc) <- as.character(sqc.header)
library(dplyr)
fam <- fam %>% select(V1, V2)
colnames(fam) <- c("FID", "IID")
fam$INDEX <- 1:nrow(fam)
sqc <- cbind(fam, sqc)
sqc %>% filter(used_in_pca_calculation == 1) %>% select(FID, IID) -> results
sqc %>% filter(used_in_pca_calculation != 1) %>% select(FID, IID) -> notresults

write.table(results, "used_in_pca_calculation.phe", quote=F, sep="\t", row.names=F, col.names=T)
write.table(notresults, "not_used_in_pca_calculation.phe", quote=F, sep="\t", row.names=F, col.names=T)
