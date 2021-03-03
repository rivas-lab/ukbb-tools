import os
import pandas as pd
import subprocess

clusters = pd.read_table('biomarkers_clusters.tsv')

qc = pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/exome/gwas/current/gwas.qc-SE02.tsv')
qc = qc[qc['#GBE_ID'].isin(list(clusters['GBE_ID']))]
qc = qc[(qc['lgc.common'] <= 3) & (qc['lgc.common'].notna()) & (qc['n_non_NA_lines'] >= 100000) & (qc['N'] >= 1000) & (qc['population'] != "e_asian") & (qc['population'] != 'metal')]

input_folder = '/oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_multitrait_exome/input/'
for cluster in clusters['cluster'].unique():
    print(cluster)
    cluster_df = clusters[clusters['cluster'] == cluster]
    paths, studies, phenos, R_phen = [], [], [], []
    for i, row in cluster_df.iterrows():
        pheno = row['GBE_ID']
        # pops = list(qc[qc['#GBE_ID'] == pheno]['population'])
        pops = ['white_british']
        for pop in pops:
            stream = os.popen('find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas -name "*.' + pheno + '.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep ' + pop)
            path = stream.read().strip()
            paths.append(path), studies.append(pop), phenos.append(pheno)
            if pop == 'white_british':
                R_phen.append('TRUE')
            else:
                R_phen.append('FALSE')
    tmp_df = pd.DataFrame(data={'path': paths, 'study': studies, 'pheno': phenos, 'R_phen': R_phen})
    tmp_df.to_csv(input_folder + str(cluster) + '.txt', sep='\t', index=False)
    stream = os.popen('sbatch -t 07-00:00:00 -p mrivas --mem=512000 --nodes=1 --cores=8 --output=mrp_logs/mrp_multitrait_exome.%j.out --wrap="/share/software/user/open/python/3.6.1/bin/python3 mrp_production.py --file ' + input_folder + str(cluster) + '.txt --R_var independent similar --variants ptv pav --sigma_m_types sigma_m_mpc_pli --filter_ld_indep --se_thresh 100 --maf_thresh 0.01 0.0005 --metadata_path /oak/stanford/groups/mrivas/ukbb24983/exome/pgen/oqfe_2020/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --out_folder /oak/stanford/groups/mrivas/users/guhan/sandbox/mrp_multitrait_exome/output/ --out_filename ' + str(cluster) + '"')
