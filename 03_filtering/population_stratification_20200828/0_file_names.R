suppressWarnings(suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
}))

get_file_names <- function(){

    file_names <- list(
        # output
        out_figs_prefix    = 'figs/sample_qc_v20200828',
        master_sqc_file    = '/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_master_sqc.20200828.phe',
        GWAS_covar_file    = '/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/ukb24983_GWAS_covar.20200828.phe',
        # input
        sqc_file           = '/oak/stanford/groups/mrivas/ukbb24983/sqc/download/ukb_sqc_v2.txt',
        sqc_colnames       = '/oak/stanford/groups/mrivas/ukbb24983/sqc/download/ukb_sqc_v2.fields.txt',
        fam_array          = '/oak/stanford/groups/mrivas/ukbb24983/fam/ukb2498_cal_v2_s488370.fam',
#         fam_exome          = '/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb_2020/ukb_exm_spb_2020.fam',
        fam_exome          = '/scratch/groups/mrivas/ukbb24983/exome-oqfe2020/ukb23155_c1_b0_v1_s200632.fam',
        remove_csv_file    = '/oak/stanford/groups/mrivas/ukbb24983/sqc/w24983_20200820.csv',
        coding1001_tsv     = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/self_reported_ethnicity/misc/coding1001.tsv',
        extracted_phe_file = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/self_reported_ethnicity/phe/ukb2007183_ukb40831_f21000.phe',
        extracted_tsv_file = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/self_reported_ethnicity/misc/ukb2007183_ukb40831_f21000.tsv',
        pop_refinement_pca = '/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/pca1_refinement',
        pop_specific_pca   = '/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification_w24983_20200828/pca2_covars',
        # semi-related
        semi_related_file  = '/oak/stanford/groups/mrivas/ukbb24983/sqc/relatedness_20200514/semi_related.fam',
        # the source of covariate
        covar_yob_tab      = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/2007183/40831/download/ukb40831.tab',
        covar_bmi_f        = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/24611/phe/INI21001.phe',
        covar_age_f        = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/time_phenotypes/misc/age_assess.phe',
        covar_CNV_f        = '/oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar_withCNV.phe',
        covar_split_dir    = '/oak/stanford/groups/mrivas/projects/degas-risk/population-split'
    )

    return(file_names)

}
