library(tidyverse)
library(data.table)
library(gridExtra)


read_great_onto_info <- function(
    assembly, assembly_version, onto, onto_version,
    onto_dir = file.path('/', 'cluster', 'data', 'ontologies')
){
    onto_term_dir <- file.path(
        onto_dir, assembly, assembly_version, onto, onto_version
    )
    
    onto_term_df <- fread(
        paste0('cat ', file.path(onto_term_dir, 'ontoTerms.canon'), '| cut -f1,2,4'), 
        data.table=FALSE
    ) %>% rename(
        Term_ID = V1,
        Desc = V2,
        DAGLvl = V3
    )
    onto_gene_df <- fread(
        file.path(onto_term_dir, 'ontoToGene.canon'), 
        data.table=FALSE
    ) %>% rename(
        Term_ID = V1,
        Gene_ID = V2
    )
    onto_info <- onto_term_df %>% inner_join(
        onto_gene_df %>% group_by(Term_ID) %>% summarise(
            n_Genes = n()
        ),
        by='Term_ID'
    )
    return(onto_info)
}

read_GREAT_res <- function(file){
    df <- fread(
        paste0('zcat ', file, ' | cut -f1,2,6,8,9', ' | sed -e "s/^#Phe/GBE_ID/g"'), 
        data.table=FALSE
    ) %>%
    rename(
        Term_ID = ID,
    ) %>%
    mutate(
        log10BPval = -log10(BPval)
    )
    return(df)
}

read_filtered_tbl <- function(
    res_dir,
    assembly, assembly_version, onto, onto_version,
    onto_n_min, onto_n_max, BFold_min, top_n, trunc_len
){
    phe_info <- fread(
        file.path(res_dir, 'jobs.tsv'),
        col.names=c('GBE_ID', 'Phe_name')
    )
    onto_info <- read_great_onto_info(
        assembly, assembly_version, onto, onto_version
    )  
    df <- read_GREAT_res(
        file.path(res_dir, 'great', paste0(onto, '.tsv.gz'))
    ) %>%
    inner_join(
        phe_info, on='GBE_ID'
    ) %>% 
    inner_join(
        onto_info, on='Term_ID'
    ) %>%
    mutate(
        Desc_trunc = str_trunc(Desc, trunc_len, "right")
    ) %>%
    mutate(
        plot_txt = paste0(Desc_trunc, ' (', Term_ID, ')')    
    ) %>% 
    filter(onto_n_min <= n_Genes) %>%
    filter(n_Genes <= onto_n_max) %>%
    filter(BFold_min  <= BFold) %>%
    group_by(GBE_ID) %>%
    arrange(-log10BPval, -BFold) %>% 
    mutate(
        Rank = rank(-log10BPval)
    ) %>% 
    filter(Rank <= top_n) %>% 
    mutate(Ontology = onto) %>% 
    ungroup() %>% 
    arrange(GBE_ID, Rank) %>%
    select(
        Ontology, GBE_ID, Phe_name, Rank, Term_ID, Desc, BPval, BFold, n_Genes, DAGLvl, plot_txt, log10BPval
    )
    return(df)
}

plot_enrichment <- function(df, gbe_id, ontology){
    df_plot <- df %>%
    filter(GBE_ID == gbe_id) %>%
    filter(Ontology == ontology)
    
    title_str <- paste0(
        (df_plot %>% select(Phe_name) %>% first())[[1]],
        ' (',
        gbe_id,
        ', ',
        ontology,
        ')'
    )
    
    p <- df_plot %>%
    ggplot(
        aes(x = reorder(as.factor(plot_txt), -Rank), y=log10BPval, fill=BFold)
    ) + 
    geom_bar(stat = 'identity') +
    coord_flip() + 
    labs(
        title = title_str,
#         x = paste0('The top ', top_n, ' significantly enriched terms \n(', ontology, ')'),
        x = paste0('The top ', dim(df_plot)[[1]], ' terms'),
        y = '-log10(p-value) for GREAT Binomial test',
        fill = 'GREAT Binomial fold'
    ) + 
    theme(legend.position="bottom")
    return(p)
}

