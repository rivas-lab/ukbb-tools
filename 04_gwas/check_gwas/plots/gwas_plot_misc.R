############################################
# GWAS plots helper functions
#
#   Yosuke Tanigawa
#   2020/9/23
############################################

suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(latex2exp)
}))


read_annotation_tbl <- function(annot_f){
    annot.arr <- fread(annot_f, sep='\t', data.table=FALSE) %>%
    mutate(
        MAF=pmin(freq, 1-freq),
        variant = paste(CHROM, POS, REF, ALT, sep=':'),
        is_outside_of_MHC = (as.numeric(CHROM) == 6 & as.numeric(POS) < 25477797) | ( as.numeric(CHROM) == 6 & 36448354 < as.numeric(POS)) | as.numeric(CHROM) != 6
    )
    PTVs <- c("frameshift_variant","splice_donor_variant","stop_gained","stop_lost","start_lost","splice_acceptor_variant")
    PAVs <- c("splice_region_variant","missense_variant","inframe_insertion","inframe_deletion")
    
    annot.arr$Csq[!(annot.arr$Consequence %in% c(PTVs,PAVs))] = "non-coding"
    annot.arr$Csq[annot.arr$Consequence %in% PAVs] = "protein-altering"
    annot.arr$Csq[annot.arr$Consequence %in% PTVs] = "protein-truncating"
    annot.arr 
}

qq_plot <- function(gwas_df){
    # generate qq plot.
    # the input data frame (gwas_df) need to have P-value column, P
    qq_p_obs <- gwas_df %>% mutate(log10P = -log10(as.numeric(P))) %>%
    arrange(-log10P) %>% drop_na(log10P) %>% pull(log10P)
    qq_p_exp <- -(qq_p_obs %>% length() %>% ppoints() %>% log10())    
    
    data.frame(Observed = qq_p_obs, Expected = qq_p_exp) %>%
    ggplot(aes(x = Expected, y = Observed)) + geom_point() + 
    geom_abline(slope=1, intercept=0, color='red') + 
    theme_bw() +  labs(
        y = latex2exp::TeX('Observed $-\\log_{10}(P)$'),
        x = latex2exp::TeX('Expected $-\\log_{10}(P)$')
    )
}

get_chrom_order <- function(){
    data.frame(CHROM = c(1:22, 'X', 'Y', 'MT'), CHROM_order=1:25, stringsAsFactors=F)
}

compute_gwas_plot_df <- function(gwas_df){
    gwas_df %>% group_by(CHROM) %>%
    summarise(chr_len = max(POS), .groups = 'drop') %>%
    left_join(get_chrom_order(), by='CHROM') %>%
    arrange(CHROM_order) %>%
    mutate(tot= cumsum(as.numeric(chr_len)) - chr_len) %>% 
    left_join(gwas_df, by="CHROM") %>%
    mutate(POScum = POS + tot) %>%
    select(-chr_len, -tot, -CHROM_order)
}

compute_x_axis_df <- function(don){
    don %>% group_by(CHROM) %>% 
    summarize(center=( max(POScum) + min(POScum) ) / 2, .groups = 'drop') %>%
    mutate(CHROM_plot = if_else(CHROM %in% c(1:11, 13, 15, 18, 21, 'X', 'Y'), as.character(CHROM), ''))
}

plot_manhattan <- function(don, pval_thr=5e-8){
    don %>% compute_x_axis_df() -> axisdf 
    don %>% ggplot( aes(x=POScum, y=-log10(P), label=repel_label) ) +
    geom_point( aes(color=as.factor(color)), alpha=0.8, size=1.3) +
    geom_hline(yintercept=-log10(pval_thr), color='red', linetype="dashed") + 
    ggrepel::geom_text_repel(size=3) +
    scale_x_continuous(label = axisdf$CHROM_plot, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5) ) +
    theme_bw() + theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + labs(
        x = 'Genomic position (chromosome)',
        y = latex2exp::TeX('$-\\log_{10\\,}P$')
    )    
}

get_Csq_color <- function(){
    data.frame(
        Csq = c('ptv', 'pav', 'pcv', 'utr', 'intron', 'others'),
        Csq_str = c('Protein-truncating variant', 'Protein-altering variant', 'Protein-coding variant', 'UTR-region variant', 'Intronic variant', 'Others'),
        Csq_color = c('#D55E00', '#E69F00', '#56B4E9', '#009E73', '#F0E442', "#A0A0A0"),
        stringsAsFactors=F
    )
}

plot_lake <- function(don){
    Csq_color <- get_Csq_color()
    don %>% compute_x_axis_df() -> axisdf 
    don %>%
    left_join(Csq_color, by='Csq') %>% select(-Csq) %>% rename('Csq'='Csq_str') %>%
    ggplot( aes(x=POScum, y=BETA, label=repel_label) ) +
    geom_point( aes(color=as.factor(Csq)), alpha=0.8, size=1.3) +
    geom_hline(yintercept=0, color="#A0A0A0") + 
    ggrepel::geom_text_repel(size=3) +
    scale_x_continuous(label = axisdf$CHROM_plot, breaks= axisdf$center) +
    theme_bw() + theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + labs(
        color = 'Consequence',
        x = 'Genomic position (chromosome)',
        y = 'BETA'
    ) +
    scale_color_manual(
        breaks = Csq_color %>% pull(Csq_str), 
        values = Csq_color %>% pull(Csq_color)
    )
}

