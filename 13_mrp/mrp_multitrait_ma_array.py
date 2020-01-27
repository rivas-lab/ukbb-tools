import os
import pandas as pd
import subprocess

clusters_to_run = pd.read_table('clusters_to_run.tsv')
clusters = pd.read_table('clusters.tsv')
clusters = clusters[clusters['cluster'].isin(list(clusters_to_run['cluster']))]
clusters = clusters.groupby('cluster')['PHEN'].apply(list).reset_index(name='phen_cluster')

sumstat_paths = pd.read_table('sumstat_paths.tsv')

for i, row in clusters.iterrows():
    paths = [], studies = [], phenos = [], R_phen = []
    phenos_in_cluster = row['phen_cluster']
    for pheno in phenos_in_cluster:
        for study in ['white_british', 'non_british_white', 'african', 'e_asian', 's_asian']:
            stream = os.popen('find /oak/stanford/groups/mrivas/ukbb24983/cal/gwas -name "*.$GBE_ID.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep ' + study)
            path = stream.read()
            paths.append(path)
            studies.append(study)
            phenos.append(pheno)
            if study == 'white_british':
                R_phen.append("TRUE")
            else:
                R_phen.append("FALSE")
    tmp_df = pd.DataFrame(data={'path': paths, 'study': studies, 'pheno': phenos, 'R_phen': R_phen})
    tmp_df.to_csv('tmp.txt', sep='\t', index=False)
    print(' '.join(row['phen_cluster']))
    stream = os.popen('sbatch -t 04:00:00 -p normal,owners,mrivas --wrap="/share/software/user/open/python/2.7.13/bin/python mrp_production.py --file tmp.txt --R_var independent similar --variants ptv pav --metadata_path /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep.tsv --out_folder /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_multitrait_array/"')
