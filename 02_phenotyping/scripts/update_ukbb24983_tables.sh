#!/bin/bash

rclone copy gdrive:rivas-lab/ukbb24983/UKBB24983_Tables.xlsx . && mv UKBB24983_Tables.xlsx ../tables
Rscript update_ukbb24983_tables.R
