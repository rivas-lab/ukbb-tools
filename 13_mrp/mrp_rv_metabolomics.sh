#!/bin/bash

# load sofware, dump which versions are used
export MODULEPATH="/home/groups/mrivas/.modules:${MODULEPATH}"
ml load python/3.6.1

phe=$1

output_folder="/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_metabolomics/"
echo -e "path\tstudy\tpheno\tR_phen" > $output_folder/$phe.tmp.txt

for POP in african e_asian s_asian white_british non_british_white related others; do
    FILEPATH=$(find /oak/stanford/groups/mrivas/users/ytanigaw/data/.20200605_NGH_UKB_phase1_metabolomics/gwas_exome/20201104/$POP -name "*.$phe.*gz" | grep -v freeze | grep -v old | grep -v ldsc);
    echo -e "$FILEPATH\t$POP\t$phe\tTRUE" >> $output_folder/$phe.tmp.txt;
done

cat $output_folder/$phe.tmp.txt

/share/software/user/open/python/3.6.1/bin/python3 mrp_production.py --file $output_folder/$phe.tmp.txt --build hg38 --R_study independent similar --R_var independent similar --variants ptv pav --filter_ld_indep --se_thresh 100 --maf_thresh 0.01 0.0005  --metadata_path /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --out_folder $output_folder

rm $output_folder/$phe.tmp.txt
