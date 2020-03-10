## Create a file for each first occurrence phenotype from UK Biobank containing y (age at event/censoring/death), status of disease (0/1) and inc (prevalent = 0, incident case = 1)
suppressMessages(require("tidyverse"))
suppressMessages(require("data.table"))
suppressMessages(require("dplyr"))
suppressMessages(require("tidyr"))
suppressMessages(require("lubridate"))


# File that has FID (f.eid) and dates of first occurrences of all phenotypes
print("Reading in first occurrences file...")
event_f <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2005693/37855/download/ukb37855.tab'
event_df <- fread(
    cmd = paste0('cat ', event_f, ' | cut -f1,', paste0(seq(14960, 17213, by=2), collapse=",")), sep='\t', data.table=T
) %>% rename("FID" = "f.eid")

# FID and dob
print("Reading in date of birth...")
dob_f <- "/oak/stanford/groups/mrivas/users/justesen/projects/disease_progression/phe/phenotype_file.phe"
dob_df <- fread(
    cmd = paste0('cat ', dob_f, ' | cut -f1,3,4'), sep='\t',data.table=T
)

# Date of censoring summer 2019? http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=42000
# Calculate age at censoring (2019-05-01)
print("Calculating age at censoring...")
dob_df$age_censo <- time_length(difftime(ISOdate(2019, 5, 1), dob_df$dob), "years")

# Age of death
print("Reading in age of death file")
death_f <- "/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/39887/download/ukb39887.tab"
death_df <- fread(
    cmd = paste0('cat ', death_f, ' | cut -f1,5914'), sep='\t', data.table=T
) %>% rename(
    FID = f.eid,
    age_death = f.40007.0.0
) 

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
        
        field_df$event_date <- as.Date(field_df[[field]], format = "%Y-%m-%d") 

        # Dates to change
        #1902-02-02	Code has event date matching participant's date of birth
        #1903-03-03	Code has event date after participant's date of birth and falls in the same calendar year as date of birth
        #2037-07-07	Code has event date in the future and is presumed to be a place-holder or other system default
        field_df <- field_df %>% 
            mutate(event_date = na_if(event_date, "2037-07-07")) 
             
        setDT(field_df)           ## converts to data.table by reference, no need for `<-`
        
        # Change 1902-02-02 to date of birth
        field_df[ event_date == "1902-02-02", event_date := as.Date(dob) ]
        # Change 1903-03-03 to date of birth + 6 months
        field_df[ event_date == "1903-03-03", event_date := as.Date(dob)%m+% months(6) ]
        
        # age at event calculation
        field_df$age_event <- time_length(difftime(as.Date(field_df$event_date), as.Date(field_df$dob)), unit = "years")

        # Set age event at birth to one month old (snpnet-cox cannot calculate based on age 0)
        field_df$age_event <- ifelse(field_df$age_event <= 0, 0.08, field_df$age_event)
        
        # Make status column
        field_df$coxnet_status_ <- ifelse(is.na(field_df$event_date), 0, NA)
        field_df$coxnet_status_ <- ifelse(!is.na(field_df$event_date), 1, field_df$coxnet_status_)

        ## Make y (age at event or age at censoring/death) column
        field_df$coxnet_y_ <- ifelse(is.na(field_df$age_event), field_df$age_censo, NA)
        field_df$coxnet_y_ <- ifelse(!is.na(field_df$age_death), field_df$age_death, field_df$coxnet_y_)
        field_df$coxnet_y_ <- ifelse(!is.na(field_df$age_event), field_df$age_event, field_df$coxnet_y_)
        
        #Seperate prevalent and incident cases (from ukbb baseline at first assessment = age0)
        field_df$age_diff <- field_df$age_event - field_df$age0

        # Variable inc for prevalent cases =0 and incident cases = 1 
        field_df$coxnet_inc_ <- ifelse(field_df$age_diff < 0, 0, NA)
        field_df$coxnet_inc_ <- ifelse(field_df$age_diff >= 0, 1, field_df$coxnet_inc_)

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