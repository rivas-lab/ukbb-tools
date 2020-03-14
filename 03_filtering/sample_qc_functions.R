library(tidyverse)
library(data.table)

read_sqc <- function(sqc_file, sqc_colnames){
    # This function reads sample QC (sqc) file
    sqc_cols <- fread(sqc_colnames, header=F, data.table=F, col.names = 'col') %>% 
    filter(col != '') %>% pull()

    fread(sqc_file, header=F, data.table=F, col.names = sqc_cols)
}

read_fam <- function(fam_file){
    # This function reads PLINK fam file
    fam_file %>% fread(
        header=F, data.table=F, 
        col.names = c('FID', 'IID', 'father', 'mother', 'sex_code', 'batch')
    ) %>%
    left_join(
        data.frame(
            sex_code = c(1, 2, 0),
            sex = c('male', 'female', 'unknown')
        ),
        by='sex_code'
    )    
}

read_remove <- function(remove_csv_file){
    # read list of individuals that should be removed from the analysis
    fread(remove_csv_file, header=F, data.table=F, col.names=c('FID')) %>%
    mutate(IID=FID, in_remove_file=TRUE)    
}

read_coding1001 <- function(coding1001_tsv){
    df <- fread(coding1001_tsv, header=T, data.table=F, sep='\t')    
    
    df %>%
    mutate(parent_id = if_else(parent_id == 0, coding, parent_id)) %>%
    left_join(
        df %>% select(coding, meaning) %>% rename(
            parent_id = coding, parent_label = meaning
        ),
        by='parent_id'
    ) %>% select(-node_id, -selectable, -parent_id) %>%
    rename(
        coding1001 = coding,
        f21000_top_label = parent_label,
        f21000_sub_label = meaning
    )    
}

read_self_reported_ethnicity <- function(extracted_phe_file){
    fread(
        cmd=paste0('cat ', extracted_phe_file, ' | sed -e "s/^#//g"'), 
        header=T, data.table=F
    )
}

read_self_reported_ethnicity_all <- function(extracted_tsv_file){
    extracted_tsv_file %>%
    fread(header=T, data.table=F) %>%
    rename(
        f21000_0 = '21000.0.0',
        f21000_1 = '21000.1.0',
        f21000_2 = '21000.2.0'    
    )
}

show_counts_for_self_reported_ethnicity <- function(master_sqc_df){
    count.wo.filter <- master_sqc_df %>% 
    count(f21000_top_label, f21000_sub_label, f21000)

    count.QC        <- master_sqc_df %>% filter(pass_QC_filter) %>% 
    count(f21000_top_label, f21000_sub_label, f21000) %>% rename(n_QC = n)

    count.QC.PCA    <- master_sqc_df %>% filter(pass_filter) %>%
    count(f21000_top_label, f21000_sub_label, f21000) %>% rename(n_QC_PCA = n)
    
    count.wo.filter %>%
    left_join(count.QC %>% select(f21000, n_QC), by='f21000')%>%
    left_join(count.QC.PCA %>% select(f21000, n_QC_PCA), by='f21000') %>%
    arrange(paste0(f21000))
}

read_master_sqc <- function(file_names){
    # read the relevant files
    sqc_df       <- read_sqc(file_names$sqc_file, file_names$sqc_colnames)
    fam_df       <- read_fam(file_names$fam_array)
    fam_exome_df <- read_fam(file_names$fam_exome)
    remove_df    <- read_remove(file_names$remove_csv_file)
    coding1001   <- read_coding1001(file_names$coding1001_tsv)
    self_reported_ethnicity     <- read_self_reported_ethnicity(file_names$extracted_phe_file)
    self_reported_ethnicity_all <- read_self_reported_ethnicity_all(file_names$extracted_tsv_file)
    
    # combine files
    bind_cols(fam_df, sqc_df) %>%
    mutate(fam_order=1:n()) %>%
    left_join(fam_exome_df %>% select(IID) %>% mutate(has_exome = T), by='IID') %>%    
    left_join(remove_df, by=c('FID', 'IID')) %>%
    left_join(self_reported_ethnicity, by=c('FID', 'IID')) %>%
    left_join(coding1001 %>% rename(f21000 = coding1001), by='f21000') %>%
    left_join(self_reported_ethnicity_all, by='IID') %>%
    replace_na(list(in_remove_file=F, has_exome = F)) %>%
    mutate(
        pass_QC_filter = (
            (putative_sex_chromosome_aneuploidy == 0) &
            (het_missing_outliers == 0) & 
            (excess_relatives == 0) & 
            (FID >= 0) &
            (IID >= 0) &
            (! in_remove_file) 
        ),
        pass_filter = (
            pass_QC_filter &
            (used_in_pca_calculation == 1)
        ),        
        self_reported_NBW = if_else(
            (is.na(f21000_0) | (f21000_0 %in% c(1, 1002, 1003))) &
            (is.na(f21000_1) | (f21000_1 %in% c(1, 1002, 1003))) &
            (is.na(f21000_2) | (f21000_2 %in% c(1, 1002, 1003))) &
            (! (is.na(f21000_0) & is.na(f21000_1) & is.na(f21000_2))),            
            T, F
        ),
        self_reported_White = if_else(
            (is.na(f21000_0) | (f21000_0 %in% c(1, 1001, 1002, 1003))) &
            (is.na(f21000_1) | (f21000_1 %in% c(1, 1001, 1002, 1003))) &
            (is.na(f21000_2) | (f21000_2 %in% c(1, 1001, 1002, 1003))) &
            (! (is.na(f21000_0) & is.na(f21000_1) & is.na(f21000_2))),            
            T, F
        ),                
        self_reported_Asian = if_else(
            (is.na(f21000_0) | (f21000_0 %in% c(3, 3001, 3002, 3003, 3004))) &
            (is.na(f21000_1) | (f21000_1 %in% c(3, 3001, 3002, 3003, 3004))) &
            (is.na(f21000_2) | (f21000_2 %in% c(3, 3001, 3002, 3003, 3004))) &
            (! (is.na(f21000_0) & is.na(f21000_1) & is.na(f21000_2))),            
            T, F
        ),                
        self_reported_Black = if_else(
            (is.na(f21000_0) | (f21000_0 %in% c(4, 4001, 4002, 4003))) &
            (is.na(f21000_1) | (f21000_1 %in% c(4, 4001, 4002, 4003))) &
            (is.na(f21000_2) | (f21000_2 %in% c(4, 4001, 4002, 4003))) &
            (! (is.na(f21000_0) & is.na(f21000_1) & is.na(f21000_2))),            
            T, F
        ),                
        self_reported_Mixed_or_Other = if_else(
            (is.na(f21000_0) | (f21000_0 %in% c(6, 2, 2001, 2002, 2003, 2004))) &
            (is.na(f21000_1) | (f21000_1 %in% c(6, 2, 2001, 2002, 2003, 2004))) &
            (is.na(f21000_2) | (f21000_2 %in% c(6, 2, 2001, 2002, 2003, 2004))) &
            (! (is.na(f21000_0) & is.na(f21000_1) & is.na(f21000_2))),            
            T, F            
        )
    ) %>%
    arrange(fam_order)
}

plot_pca_add_threshold <- function(p){
    p + 
    
    geom_vline(xintercept=260) +  # African      260 <= PC1        &&   50 <= PC2
    geom_hline(yintercept=50) +  

    geom_vline(xintercept=40) +   # South Asian   40 <= PC1 <= 120 && -170 <= PC2 <= -80
    geom_vline(xintercept=120) +    
    geom_hline(yintercept=-80) +  
    geom_hline(yintercept=-170) +

    geom_vline(xintercept=130) +   # East Asian   130 <= PC1 <= 170 &&         PC2 <= -230
    geom_vline(xintercept=170) +        
    geom_hline(yintercept=-230) + 

    geom_vline(xintercept=-20) +  # European     -20 <= PC1 <=  40 &&  -25 <= PC2 <= 10
    geom_hline(yintercept=-25) +  
    geom_hline(yintercept=10)    
}

plot_pca_self_reported <- function(df){    
    df %>% 
    drop_na(f21000) %>%    
    filter(pass_filter, f21000 != -3, f21000 != -1) %>%
    ggplot(aes(x = PC1, y = PC2, color=as.factor(f21000_top_label))) %>%
    plot_pca_add_threshold() +        
    geom_point(alpha=.025) +    
    theme_bw() + 
    theme(legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(alpha = 1), nrow=2, byrow=TRUE)) +
    labs(title='Global PCs', color='Self-reported ethnicity')    
}

plot_pca_population <- function(df, x_axis, y_axis){    
    df %>% 
    drop_na(f21000) %>% 
    filter(pass_filter, f21000 != -3, f21000 != -1) %>%
    mutate(
        population=if_else(
            (
                (population == 's_asian_outlier') | 
                (population == 'e_asian_outlier')
            ), '', population
        ),
        population=na_if(population, '')
    ) %>%   
    replace_na(list(population='others')) %>%
    rename(plot_x = x_axis, plot_y = y_axis) %>%
    ggplot(aes(x = plot_x, y = plot_y, color=as.factor(f21000_top_label))) +
    geom_point(alpha=.025) +
    theme_bw() + 
    theme(legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(alpha = 1), nrow=2, byrow=TRUE)) +
    labs(x = x_axis, y = y_axis, title='Global PCs', color='Self-reported ethnicity') + 
    facet_wrap( ~ population, ncol=3) 
}


define_populations <- function(sqc_df){
    pop_df <- sqc_df %>% filter(pass_filter) %>% 
    mutate(
        l_african  = if_else(
            260 <= PC1 &                50 <= PC2 &
            (! self_reported_Asian) & (! self_reported_White) & (! self_reported_Mixed_or_Other), T, F),
        l_s_asian  = if_else(
             40 <= PC1 & PC1 <= 120 & -170 <= PC2 & PC2 <= -80 &
            (! self_reported_Black) & (! self_reported_White) & (! self_reported_Mixed_or_Other), T, F),
        l_e_asian  = if_else(
            130 <= PC1 & PC1 <= 170 &               PC2 <= -230 &
            (! self_reported_Black) & (! self_reported_White) & (! self_reported_Mixed_or_Other), T, F),
        l_European = if_else(
            -20 <= PC1 & PC1 <=  40 &  -25 <= PC2 & PC2 <= 10, T, F),
        l_non_british_white = if_else(l_European & self_reported_NBW, T, F),
        l_white_british  = if_else(l_European & in_white_British_ancestry_subset, T, F)
    ) %>% 
    select(IID, l_white_british, l_non_british_white, l_african, l_s_asian, l_e_asian) %>% 
    gather("population", "label", -IID) %>% filter(label) %>% select(-label) %>%
    mutate(population = str_replace(population, '^l_', ''))
    
    sqc_df %>% left_join(pop_df, by='IID') %>%
    arrange(fam_order)
}

show_population_counts <- function(df){
    df %>% 
    drop_na(population) %>%
    count(population, genotyping_array) %>% 
    spread(genotyping_array, n, fill= 0) %>%
    left_join(
        df %>% 
        drop_na(population) %>%
        count(population, has_exome) %>% 
        mutate(has_exome = paste0(has_exome)) %>% 
        spread(has_exome, n, fill= 0) %>%
        rename(w_exome = 'TRUE', wo_exome = 'FALSE') %>%
        mutate(n = w_exome + wo_exome) %>%
        select(population, w_exome, wo_exome, n),    
        by='population'
    ) %>%
    arrange(-n)    
}


read_eigenvec <- function(dir_name, pops){
    eigenvec <- list()
    for (pop in pops){
        if(pop != 'white_british'){
            eigenvec[[pop]] <- file.path(dir_name, paste0('ukb24983_', pop, '_pca.eigenvec')) %>%
            fread(header=T, data.table=F) %>% rename(FID = '#FID') %>% mutate(population = pop)  
        }
    }
    eigenvec %>% bind_rows()
}

apply_threshold <- function(eigenvec_df, x_axis, y_axis, pop, x_lim, y_lim){

    x_min <- x_lim[1]
    x_max <- x_lim[2]
    y_min <- y_lim[1]
    y_max <- y_lim[2]

    p <- eigenvec_df %>% 
    filter(population == pop) %>%
    left_join(master_sqc_pop_df %>% select(IID, f21000, f21000_top_label), by='IID') %>%
    drop_na(f21000) %>% filter(f21000 != -3, f21000 != -1) %>%
    rename(plot_x = x_axis, plot_y = y_axis) %>%
    ggplot(aes(x = plot_x, y = plot_y, color=as.factor(f21000_top_label))) +
    geom_point(alpha=.025) + 
    theme_bw() + theme(legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(alpha = 1), nrow=2, byrow=TRUE)) +
    labs(
        x = x_axis, y = y_axis, 
        title='PCs within each population',
        color='Self-reported ethnicity'
    ) + 
    geom_vline(xintercept=x_min) +
    geom_vline(xintercept=x_max) +
    geom_hline(yintercept=y_min) +
    geom_hline(yintercept=y_max)

    n_before_filter <- eigenvec_df %>% 
    filter(population == pop) %>% 
    nrow()

    n_after_filter <- eigenvec_df %>% 
    rename(plot_x = x_axis, plot_y = y_axis) %>%
    filter(
        population == pop,
        x_min <= plot_x, plot_x <= x_max,
        y_min <= plot_y, plot_y <= y_max    
    ) %>% 
    nrow()
    
    print(sprintf(
        'Number of individuals: %d (before filter) --> %d (after filter)',
        n_before_filter, n_after_filter
    ))
    
    p
}

population_def_refinement <- function(master_sqc_pop_df, eigenvec_df){
    population_refinement_df <- eigenvec_df %>% 
    mutate(
        population_refined_s_asian = if_else(
            (population == 's_asian') &
            (-0.02 <= PC1) & (PC1 <= 0.03) &
            (-0.05 <= PC2) & (PC2 <= 0.02), T, F
        ),
        population_refined_e_asian = if_else(
            (population == 'e_asian') &
            (-0.01 <= PC1) & (PC1 <= 0.02) &
            (-0.02 <= PC2) & (PC2 <= 0), T, F
        ),    
    ) %>% 
    select(IID, population_refined_s_asian, population_refined_e_asian)

    WB_IIDs <- master_sqc_pop_df %>% 
    filter(population == 'white_british') %>%
    select(IID) %>% pull()

    master_sqc_pop_df %>% 
    mutate(
        population = if_else(population == 's_asian', 's_asian_outlier', population),
        population = if_else(population == 'e_asian', 'e_asian_outlier', population)
    ) %>% 
    left_join(population_refinement_df, by='IID') %>%
    mutate(
        population = if_else(population_refined_s_asian, 's_asian', population),
        population = if_else(population_refined_e_asian, 'e_asian', population),
        population = if_else((FID %in% WB_IIDs), 'white_british', population)
    ) %>%
    select(-population_refined_s_asian, -population_refined_e_asian) %>% 
    arrange(fam_order)
}

plot_local_pc <- function(evec_df, master_sqc_pop_ref_df, x_axis, y_axis){
    evec_df %>% 
    select(-population) %>%
    inner_join(
        master_sqc_pop_ref_df %>% 
        mutate(
            population=if_else(
                (
                    (population == 's_asian_outlier') | 
                    (population == 'e_asian_outlier')
                ), '', population
            ),
            population=na_if(population, '')
        ) %>%
        drop_na(population, f21000) %>%
        filter(f21000 != -3, f21000 != -1) %>%
        select(IID, f21000, f21000_top_label, population), 
        by='IID'
    ) %>%
    rename(plot_x = x_axis, plot_y = y_axis) %>%
    ggplot(aes(x = plot_x, y = plot_y, color=as.factor(f21000_top_label))) +
    geom_point(alpha=.025) + 
    theme_bw() + 
    theme(legend.position="bottom") +
    guides(colour = guide_legend(override.aes = list(alpha = 1), nrow=2, byrow=TRUE)) +
    labs(
        x = x_axis, y = y_axis, 
        title='PCs within each population',
        color='Self-reported ethnicity'
    ) +
    facet_wrap( ~ population, ncol=2)     
}

replace_colnames <- function(df, pattern, replacement){
    new_colnames <- df %>% 
    colnames() %>%
    lapply(
        function(x){
            str_replace(x, pattern, replacement)
        }
    )    
    colnames(df) <- new_colnames
    df
}
