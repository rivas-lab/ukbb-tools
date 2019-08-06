sqc.header <- read.table("/oak/stanford/groups/mrivas/ukbb24983/sqc/download/ukb_sqc_v2.fields.txt", header=F, stringsAsFactors=F)$V1
fam <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2.fam", header=F)
sqc <- data.table::fread("/oak/stanford/groups/mrivas/ukbb24983/sqc/download/ukb_sqc_v2.txt", header=F)
colnames(sqc) <- as.character(sqc.header)
library(dplyr)
fam <- fam %>% select(V1, V2)
colnames(fam) <- c("FID", "IID")
fam$INDEX <- 1:nrow(fam)
sqc <- cbind(fam, sqc)

sqc %>% filter(dQC > 0.9, Internal_Pico_ng_uL > 0.1, Submitted_Gender == Inferred_Gender, het_missing_outliers == 0, putative_sex_chromosome_aneuploidy == 0, excluded_from_kinship_inference == 0) %>% select(FID, IID) -> results
#sqc %>% filter(dQC > 0.9, Internal_Pico_ng_uL > 0.1, Submitted_Gender == Inferred_Gender, het_missing_outliers == 0, putative_sex_chromosome_aneuploidy == 0, excluded_from_kinship_inference == 0, heterozygosity_pc_corrected < 0.05) %>% select(FID, IID) -> results

write.table(results, "pass_qc.phe", quote=F, sep="\t", row.names=F, col.names=T)

