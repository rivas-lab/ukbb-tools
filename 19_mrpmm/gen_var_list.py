#!/bin/python

import pandas as pd

metadata = pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv')

#genes_of_interest = ['PCSK9', 'LMF1', 'LIPC', 'LDLRAP1', 'LDLR', 'CETP', 'APOC3', 'APOB', 'ABCG8', 'ABCG5', 'ABCA1']
creat = pd.read_table('creat_genes_to_run.tsv')
genes_of_interest = list(creat['gene'])

ptv =  [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "start_lost",
        "stop_lost",
]


metadata = metadata[(metadata['gene_symbol'].isin(genes_of_interest)) & (metadata['ld_indep'] == True) & (metadata['most_severe_consequence'].isin(ptv))]
metadata = metadata[['V']]
#metadata.to_csv('HDL_LDL_TG_variants.lst', sep='\t', index=False)
metadata.to_csv('creat_variants.lst', sep='\t', index=False)
