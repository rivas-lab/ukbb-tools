import os
import pandas as pd
import subprocess

clusters_to_run = pd.read_table('clusters_to_run.tsv')
clusters = pd.read_table('clusters.tsv')
clusters = clusters[clusters['cluster'].isin(list(clusters_to_run['cluster']))]
clusters.to_csv('clusters_subset.tsv', sep='\t', index=False)
clusters = clusters.groupby('cluster')['PHEN'].apply(list).reset_index(name='phen_cluster')

sumstat_paths = pd.read_table('sumstat_paths.tsv')
input_folder = '/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_multitrait_array/input/'

for i, row in clusters.iterrows():
    paths, phenos = [], []
    phenos_in_cluster = row['phen_cluster']
    for pheno in phenos_in_cluster:
        stream = os.popen('find /oak/stanford/groups/mrivas/ukbb24983/cal/gwas -name "*.' + pheno + '.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep white_british')
        path = stream.read().strip()
        if path != "":
            paths.append(path)
            phenos.append(pheno)
        else:
            print("Summary statistic file for " + pheno + "not found!")
    tmp_df = pd.DataFrame(data={'path': paths, 'study': ['white_british'] * len(paths), 'pheno': phenos, 'R_phen': ['TRUE'] * len(paths)})
    tmp_df.to_csv(input_folder + str(row['cluster']) + '.txt', sep='\t', index=False)
    print(' '.join(row['phen_cluster']))
    stream = os.popen('sbatch -t 04:00:00 -p normal,owners,mrivas --mem=32000 --wrap="/share/software/user/open/python/2.7.13/bin/python mrp_production.py --file ' + input_folder + str(row['cluster']) + '.txt --R_var independent similar --variants ptv pav --metadata_path /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --out_folder /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_multitrait_array/output/ --out_filename ' + str(row['cluster']) + '"')
