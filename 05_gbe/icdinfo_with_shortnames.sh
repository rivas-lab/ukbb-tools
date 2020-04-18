#!/bin/bash

rclone copy gdrive:rivas-lab/ukbb24983/GBE_names.xlsx .
for pop in white_british non_british_white s_asian e_asian african; do bash update_icdinfo.sh ${pop}; done
Rscript icdinfo_with_shortnames.R
for pop in non_british_white s_asian e_asian african; do rm icdinfo.${pop}.txt; done
rclone copy GBE_names.xlsx gdrive:rivas-lab/ukbb24983
