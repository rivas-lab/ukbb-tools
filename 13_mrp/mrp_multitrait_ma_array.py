import os
import pandas as pd
import subprocess

clusters_to_run = [85]#pd.read_table('clusters_to_run.tsv')
clusters = pd.read_table('clusters.tsv')
clusters = clusters[clusters['cluster'].isin(clusters_to_run)] #(list(clusters_to_run['cluster']))]
#clusters.to_csv('clusters_subset.tsv', sep='\t', index=False)
clusters = clusters.groupby('cluster')['PHEN'].apply(list).reset_index(name='phen_cluster')

n_thresh = 100

sumstat_paths = pd.read_table('sumstat_paths.tsv')
input_folder = '/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_multitrait_ma_array/input/'
pops = ['white_british', 'non_british_white', 'african', 'e_asian', 's_asian']
pop_cols = [7, 8, 9, 10, 11]
phen_info = pd.read_table('../05_gbe/phenotype_info.tsv')

for i, row in clusters.iterrows():
    paths, studies, phenos, R_phen = [], [], [], []
    phenos_in_cluster = row['phen_cluster']
    for pheno in phenos_in_cluster:
        subset_phen = phen_info[phen_info['#GBE_ID'] == pheno]
        pop_nums = subset_phen[subset_phen.columns[pop_cols]].values.tolist()[0]
        for pop, pop_num in zip(pops, pop_nums):
            if pop_num >= n_thresh:
                stream = os.popen('find /oak/stanford/groups/mrivas/ukbb24983/array-combined/gwas -name "*.' + pheno + '.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep ' + pop)
                path = stream.read().strip()
                if path != "":
                    paths.append(path)
                    studies.append(pop)
                    phenos.append(pheno)
                    if pop == 'white_british':
                        R_phen.append('TRUE')
                    else:
                        R_phen.append('FALSE')
                else:
                    print("Summary statistic file for " + pheno + "not found!")
            else:
                print(pheno + " " + pop + " does not meet minimum N threshold!")
    tmp_df = pd.DataFrame(data={'path': paths, 'study': studies, 'pheno': phenos, 'R_phen': R_phen})
    tmp_df.to_csv(input_folder + str(row['cluster']) + '.txt', sep='\t', index=False)
    print(' '.join(row['phen_cluster']))
    stream = os.popen('sbatch -t 04:00:00 -p normal,owners,mrivas --mem=32000 --output=mrp_logs/mrp_multitrait_ma_array.%j.out --wrap="/share/software/user/open/python/2.7.13/bin/python mrp_production.py --file ' + input_folder + str(row['cluster']) + '.txt --R_study independent similar --R_var independent similar --variants ptv pav --metadata_path /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --out_folder /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_multitrait_ma_array/output/ --out_filename ' + str(row['cluster']) + '"')
