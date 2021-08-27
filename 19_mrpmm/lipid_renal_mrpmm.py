import pandas as pd
import os

exome = pd.read_table('mrp200kexomeresults.tsv')

def file_find(pheno):
    stream = os.popen('find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas -name "*.' + pheno + '.*gz" | grep -v freeze | grep -v old | grep -v ldsc | grep metal | grep -v plink')
    path = stream.read().strip()
    return path

for trait_set, output_folder, fout in zip([['INI30030700', 'INI10030720', 'INI10030700'], ['INI10030760', 'INI20030780', 'INI10030870']], ['/oak/stanford/groups/mrivas/users/guhan/sandbox/mrpmm_renal', '/oak/stanford/groups/mrivas/users/guhan/sandbox/mrpmm_lipid'], ['renal', 'lipid']):
    exome_genes = list(set(exome[(exome['#GBE_ID'].isin(trait_set)) & (exome['log_10_BF_study_similar_var_similar_sigma_m_mpc_pli_ptv'] >= 5)]['gene']))
    print(fout)
    print(exome_genes)
    print(len(exome_genes))
    paths = [file_find(pheno) for pheno in trait_set]
    tmp_df = pd.DataFrame(data={'path': paths, 'study': ['metal'] * len(paths), 'pheno': trait_set, 'R_phen': ["TRUE"] * len(paths)})
    tmp_df.to_csv('_'.join(trait_set) + ".tmp.txt", sep='\t', index=False)
    metadata = pd.read_table('../13_mrp/ref/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv')
    metadata[metadata['gene_symbol'].isin(exome_genes)][['V']].to_csv('_'.join(trait_set) + "_variants.tsv", sep='\t', index=False)
    cmd = 'sbatch -t 7-00:00:00 -p mrivas --mem=200000 --nodes=1 --cores=16 --wrap="python3 mrpmm.py --variants ' + '_'.join(trait_set) + "_variants.tsv" + ' --phenotypes ' + '_'.join(trait_set) + ".tmp.txt" + ' --metadata_path ../13_mrp/ref/ukb_exm_oqfe-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv --build hg38 --out_folder ' + output_folder + ' --C 1 2 3 4 5 6 7 8 --se_thresh 100 --maf_thresh 0.01 --variant_filters ptv --fout ' + fout + '"'
    print(cmd)
    stream = os.popen(cmd)
