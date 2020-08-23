#!/bin/bash
set -beEuo pipefail

long_df_f='/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc_upset/ukb37855_ukb40831_icd.tsv.zst'
coding_f='/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ukbb-phenotyping/sr_icd_map.csv'
long_df_annot_f='/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc_upset/ukb37855_ukb40831_icd.annot.tsv.gz'

ml load R/3.6 gcc

if [ ! -f "${long_df_annot_f%.gz}" ] && [ ! -f "${long_df_annot_f}" ] ; then

    Rscript 2_annotate_extracted_tab.R ${long_df_f} ${coding_f} ${long_df_annot_f%.gz}
    echo ${long_df_annot_f%.gz}

fi

# compress and index
bgzip --compress-level 9 -@6 -f ${long_df_annot_f%.gz}
echo ${long_df_annot_f}

tabix ${long_df_annot_f} -s 1 -b 2 -e 2
echo ${long_df_annot_f}.tbi
