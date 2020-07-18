args <- commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

####################################################################
x <- args[1]
out_d <- '/scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc'
# annot_f <- '/oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.6302020.tsv.gz'
annot_f <- '/scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc/variant_filter_table.6302020.mafonly.tsv.gz'
# zcat /oak/stanford/groups/mrivas/private_data/ukbb/variant_filtering/variant_filter_table.6302020.tsv.gz | cut -f 1-4,19 | bgzip -@6 >  /scratch/groups/mrivas/ukbb24983/array-combined/gwas-qc/variant_filter_table.6302020.mafonly.tsv.gz
####################################################################

lgc<-function(x){
	x <- as.numeric(x[!is.na(x)])
	n <- length(x)
	x2obs <- qchisq(x,1,lower.tail=FALSE)
  	x2exp <- qchisq(1:n/n,1,lower.tail=FALSE)
    lambda <- median(x2obs)/median(x2exp)
    return(lambda)
}

lgc_main <- function(ss, pop, fsum, annot, out_f){
	df <- merge(ss, annot, by.x = c("#CHROM","POS","REF","ALT"), by.y = c("CHROM","POS","REF","ALT"), all.x = TRUE)
	write.table(paste(fsum, pop, "common",     lgc(df$P[df$maf > .05]), sep='\t'), quote = FALSE, col.names = FALSE, row.names = FALSE, file = out_f, append = FALSE)
	write.table(paste(fsum, pop, ".05",        lgc(df$P[df$maf <= .05]), sep='\t'), quote = FALSE, col.names = FALSE, row.names = FALSE, file = out_f, append = TRUE)
	write.table(paste(fsum, pop, ".01",        lgc(df$P[df$maf <= .01]), sep='\t'), quote = FALSE, col.names = FALSE, row.names = FALSE, file = out_f, append = TRUE)
	write.table(paste(fsum, pop, ".001-.01",   lgc(df$P[df$maf <= .01 & df$maf >= .001]), sep='\t'), quote = FALSE, col.names = FALSE, row.names = FALSE, file = out_f, append = TRUE)
	write.table(paste(fsum, pop, ".0001-.001", lgc(df$P[df$maf <= .001 & df$maf >= .0001]), sep='\t'), quote = FALSE, col.names = FALSE, row.names = FALSE, file = out_f, append = TRUE)
	write.table(paste(fsum, pop, ".0001",      lgc(df$P[df$maf <= .0001]), sep='\t'), quote = FALSE, col.names = FALSE, row.names = FALSE, file = out_f, append = TRUE)
}
####################################################################

pop   <- basename(dirname(x))
fsum  <- str_split(str_split(basename(x), 'ukb24983_v2_hg19.', simplify = TRUE)[1,2],'.array', simplify = TRUE)[1,1]
out_f <- file.path(out_d, pop, sprintf('%s.%s.qc.txt', pop, fsum))

message(sprintf('%s %s %s', pop, fsum, out_f))

if(! dir.exists(file.path(out_d, pop))){
	dir.create(file.path(out_d, pop))
}
annot <- fread(annot_f, header = TRUE, sep = "\t", data.table = FALSE)
ss <- fread(x,data.table = FALSE, sep = "\t", header = TRUE)
lgc_main(ss, pop, fsum, annot, out_f)