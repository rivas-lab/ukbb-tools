require(readxl)

data <- read_excel('../tables/UKBB24983_Tables.xlsx')
write.table(data, file='../tables/ukbb24983_tables.tsv', quote=FALSE, sep='\t', row.names = F, na="")
