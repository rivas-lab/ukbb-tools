#!/usr/bin/env Rscript
#### ARG PARSER ####

# parse command line args: this controls variations of model
args <- commandArgs(trailingOnly = TRUE)

suppressMessages(require(tidyverse))
suppressMessages(require(data.table))
suppressMessages(require(reshape))

if (length(args) <= 1) {
    stop("Need at least two phenotypes to run analysis!")
} else if (length(args) >= 2) {
    print("Preparing phenotype file")
    
    # read files 
    phes_to_combine <- paste("/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/master_phe/opcs/phenotypefiles/", args, ".phe", sep='')
    
    # Merge files
    file.list <- lapply(phes_to_combine, read.table, header=TRUE)
    merged_df <- merge_recurse(file.list, all.x = TRUE, by="FID")
                        
    # Filter on disease = 1 
    merged_df <- merged_df %>% filter(get(paste("coxnet_status_", args[1], sep = "")) == 1)

    ys = c()
    statuses = c()
    sum_status_vec = c()
   
    for (row in 1:nrow(merged_df)) {
        # If all other statuses are all 0s, then take the age, status == 0, etc and add to a new dataframe
        sum_statuses <- sum(merged_df[row, paste("coxnet_status_", args[-1], sep="")])
        sum_status_vec[row] = sum_statuses
        if (sum_statuses == 0) {
            ys[row] = merged_df[row, paste("coxnet_y_", args[2], sep="")]
            statuses[row] = merged_df[row, paste("coxnet_status_", args[2], sep="")] 
        # If at least one status is 1, look at age columns across all statuses and pick phenotype with mininum age
        } else if (sum_statuses >= 1) {
            ys[row] = min(merged_df[row, paste("coxnet_y_", args[-1], sep="")])
            statuses[row] = 1
        }
    }
    
    # Create new dataframe with FID, pooled age, pooled status, etc.
    merged_df$combined_y = ys
    merged_df$combined_status = statuses
    
    print(dim(merged_df))
    
    # Filter based on t1 age vs first opcs event
    merged_df <- merged_df %>% filter(get(paste('coxnet_y_', args[1], sep="")) < combined_y) %>% 
        select(FID, combined_y, combined_status)
   
    print("after filtering on age < opcs event")
    print(dim(merged_df))
     
    # Add out_folder later - for now write new file to working directory
    #out_folder = "//"
    # Write the results to file 
    filename <- paste(paste("snpnet_cox", paste(args, collapse='_'), sep="_"), ".phe", sep="")
    write.table(merged_df, file=paste(filename, sep=""), quote=F, row.names=F, sep='\t')
    #write.table(merged_df, file=paste(out_folder, filename, sep=""), quote=F, row.names=F, sep='\t')
}   
                      
  # disease <- args[1] # First occurrence
  # opcs <- args[-1] # all other args are opcs to be combined                    