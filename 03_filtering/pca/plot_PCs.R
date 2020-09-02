suppressWarnings(suppressPackageStartupMessages({ library(tidyverse); library(data.table) }))
args <- commandArgs(trailingOnly=TRUE)
evec_f <- args[1]
evec_p <- args[2]

x <- 'PC1'
y <- 'PC2'
if(length(args)>2){ x <- args[3] ; }
if(length(args)>3){ y <- args[4] ; }

#############
evec <- fread(evec_f, colClasses=c('#FID'='character', 'IID'='character')) %>%
rename('FID'='#FID')

p1 <- evec %>%
rename(plot_x=x, plot_y=y) %>%
ggplot(aes(x=plot_x, y=plot_y)) +
stat_density_2d(aes(fill = ..level..), geom = "polygon") +
labs(title = file.path(basename(dirname(evec_f)), basename(evec_f)), x=x, y=y) +
theme_bw()

p2 <- evec %>%
rename(plot_x=x, plot_y=y) %>%
ggplot(aes(x=plot_x, y=plot_y)) +
geom_point(alpha=.05) +
labs(title = file.path(basename(dirname(evec_f)), basename(evec_f)), x=x, y=y) +
theme_bw()

ggsave(evec_p, gridExtra::arrangeGrob(p1, p2, ncol=2), width=12, height=6)
