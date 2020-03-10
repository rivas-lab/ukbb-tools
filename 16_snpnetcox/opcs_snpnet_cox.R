# This script prepares individual files for OPCS4 codes with four columns FID, coxnet_y_[OPCS_ID], coxnet_status_[OPCS_ID], coxnet_inc_[OPCS_ID]
suppressMessages(require("tidyverse"))
suppressMessages(require("data.table"))
suppressMessages(require("dplyr"))
suppressMessages(require("tidyr"))
suppressMessages(require("lubridate"))

# File that has FID and date of OPCS4
print("Reading in opcs file...")
#event_df <- read.table('/oak/stanford/groups/mrivas/users/justesen/disease_progression/phe/opcs_phe/opcs.phe', sep='\t', head=T)
event_df <- fread("/oak/stanford/groups/mrivas/users/justesen/projects/disease_progression/phe/opcs.formatting.phe", 
                 sep = "\t", header= TRUE)
# FID and dob
print("Reading in date of birth...")
dob_f <- "/oak/stanford/groups/mrivas/users/justesen/projects/disease_progression/phe/phenotype_file.phe"
dob_df <- fread(
    cmd = paste0('cat ', dob_f, ' | cut -f1,3,4'), sep='\t',data.table=T
)

# Date of censoring 
# Calculate age at censoring (2017-04-01)
print("Calculating age at censoring...")
dob_df$age_censo <- time_length(difftime(ISOdate(2017, 4, 1), dob_df$dob), "years")

# Age of death
print("Reading in age of death file")
death_f <- "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/39887/download/ukb39887.tab"
death_df <- fread(
    cmd = paste0('cat ', death_f, ' | cut -f1,5914'), sep='\t', data.table=T
) 
colnames(death_df)[1] = 'FID'
colnames(death_df)[2] = 'age_death'

# Merge
print("Merging dataframes...")
merged_df <- Reduce(function(x, y) merge(x, y, by="FID", all=TRUE), list(event_df, dob_df, death_df))

print("Memory before loop:")
print(gc())

# Loop over all phenotypes
for (field in colnames(event_df)){

    print(field)

    if (field != "FID") {
        # Select 
        field_df <- merged_df %>% select("FID", "dob", "age0", "age_censo", "age_death", field)
        
        # Make field into the event_date
        field_df$event_date <- as.Date(field_df[[field]], format = "%Y-%m-%d")
        
        # age at event calculation
        field_df$age_event <- time_length(difftime(as.Date(field_df$event_date), as.Date(field_df$dob)), unit = "years")
        
        # Set age event at birth to one month old (snpnet-cox does not allow age 0)
        field_df$age_event <- ifelse(field_df$age_event <= 0, 0.08, field_df$age_event)
        
        # Make status column
        field_df$coxnet_status_ <- ifelse(is.na(field_df$age_event), 0, NA)
        field_df$coxnet_status_ <- ifelse(!is.na(field_df$age_event), 1, field_df$coxnet_status_)

        # Make y (age at event or age at censoring/death) column
        field_df$coxnet_y_ <- ifelse(is.na(field_df$age_event), field_df$age_censo, NA)
        field_df$coxnet_y_ <- ifelse(!is.na(field_df$age_death), field_df$age_death, field_df$coxnet_y_)
        field_df$coxnet_y_ <- ifelse(!is.na(field_df$age_event), field_df$age_event, field_df$coxnet_y_) 
        
        #Seperate prevalent and incident cases
        field_df$age_diff <- field_df$age_event - field_df$age0

        # Variable inc for prevalent cases =0 and incident cases = 1 
        field_df$coxnet_inc_ <- ifelse(field_df$age_diff < 0, 0, NA)
        field_df$coxnet_inc_ <- ifelse(field_df$age_diff >= 0, 1, field_df$coxnet_inc_)

        # Select columns for final table
        status_df <- field_df %>%
            select(
                FID,
                coxnet_y_,
                coxnet_status_,
                coxnet_inc_
            )


        # rename for individual phenotype
        cols <- c("coxnet_y_", "coxnet_status_","coxnet_inc_")
        status_df <- rename_at(status_df, cols, list(~paste0(.,field)))
        
        write.table(status_df, paste0(field, ".phe"), quote=FALSE, sep="\t", row.names=FALSE)
    }
}