library(dplyr)
library(mvtnorm)
library(MCMCpack)

set.seed(42)

covar <- read.table("/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe", header=T)
globalpcs <- covar[,c(1:2,6:45)]

prefix = "lambda10"

for (pop in c("african", "chinese", "mixedwhiteafrican", "southasian", "white", "mixedwhiteasian")) {
    putative <- read.table(sprintf("self_identified_%s.phe", pop), header=F, col.names=c("FID", "IID"))
    
    local.si <- read.table(sprintf("self_identified_%s.firstpca10.phe", pop), header=T) %>% inner_join(putative)
    
    #pdf(sprintf("aberrant.self_identified_%s.pdf", pop), useDingbats=F)
    aberrant::aberrant(local.si[,c("PC1_AVG", "PC2_AVG")], 5) -> local.si.12
    #aberrant::plot(local.si.12)

    aberrant::aberrant(local.si[,c("PC3_AVG", "PC4_AVG")], 5) -> local.si.34
    #aberrant::plot(local.si.34)

    aberrant::aberrant(local.si[,c("PC5_AVG", "PC6_AVG")], 5) -> local.si.56
    #aberrant::plot(local.si.56)
    #dev.off()
    
    intersect(intersect(local.si.12$inlier, local.si.34$inlier), local.si.56$inlier) -> all.inlier
    
    local.si[all.inlier,] %>% dplyr::select(FID) %>% mutate(IID=FID) %>% write.table(sprintf("%s/inlier.self_identified_%s.phe", prefix, pop), quote=F, sep="\t", row.names=F, col.names=T)
    
    intersect(intersect(local.si.12$outlier, local.si.34$outlier), local.si.56$outlier) -> all.outlier
    local.si[all.outlier,] %>% dplyr::select(FID) %>% mutate(IID=FID) %>% write.table(sprintf("%s/outlier.self_identified_%s.phe", prefix, pop), quote=F, sep="\t", row.names=F, col.names=T)
    save(local.si.12, local.si.34, local.si.56, local.si, all.inlier, all.outlier, file=sprintf("%s/aberrant.self_identified_%s.RData", prefix, pop))
}
