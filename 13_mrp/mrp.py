from __future__ import division
import pandas as pd
from functools import partial, reduce
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import numpy.matlib as npm
from scipy.stats import describe
from scipy.stats import beta
from collections import defaultdict
import subprocess
import math
import os
import sys

#Read in summary statistics from GBE sumstats
def read_in_summary_stats(disease_string):
    findCMD = 'find /oak/stanford/groups/mrivas/ukbb24983/exome/gwas/ -name "*' + disease_string + '.*" | grep -v "exome-spb.log" | grep -v freeze | grep -v old'
    out = subprocess.Popen(findCMD,shell=True,stdin=subprocess.PIPE, 
                        stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdout, stderr) = out.communicate()
    sumstat_file = stdout.strip().decode("utf-8") 
    print("Sumstat file:")
    print(sumstat_file)
    df = pd.read_csv(sumstat_file, sep='\t')
    df.insert(loc=0, column='V', value=df['#CHROM'].astype(str).str.cat(df['POS'].astype(str), sep=':').str.cat(df['REF'], sep=':').str.cat(df['ALT'], sep=':'))
    if sumstat_file.endswith('linear.gz'):
        df = df[['V', 'ID', 'A1', 'TEST', 'OBS_CT', 'BETA', 'SE', 'T_STAT', 'P']]
    else:
        df = df[['V', 'ID', 'A1', 'FIRTH?', 'TEST', 'OBS_CT', 'OR', 'SE', 'Z_STAT', 'P']]
    return df

#Need to have: gene, consequence, pop/global allele frequencies.
def read_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, sep='\t')
    return metadata

def set_sigmas(df):
    to_filter = ['regulatory_region_variant', 'intron_variant', 'intergenic_variant', 'downstream_gene_variant', 'mature_miRNA_variant', 'non_coding_transcript_exon_variant', 'upstream_gene_variant']
    print("Before consequence filter:")
    print(len(df))
    df = df[~df.most_severe_consequence.isin(to_filter)]
    print("After consequence filter:")
    print(len(df))
    print("Setting sigmas...")
    sigma_m_ptv = 0.2
    sigma_m_pav = 0.05
    sigma_m_pc = 0.05
    ptv = ['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'start_lost', 'stop_lost']
    pav = ['protein_altering_variant', 'inframe_deletion', 'inframe_insertion', 'splice_region_variant', 'start_retained_variant', 'stop_retained_variant', 'missense_variant']
    proximal_coding = ['synonymous_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'coding_sequence_variant', 'incomplete_terminal_codon_variant', 'TF_binding_site_variant']
    sigma_m = dict([(variant, sigma_m_ptv) for variant in ptv] + [(variant, sigma_m_pav) for variant in pav]
                   + [(variant, sigma_m_pc) for variant in proximal_coding])
    category_dict = dict([(variant, 'ptv') for variant in ptv] + [(variant, 'pav') for variant in pav]
                   + [(variant, 'proximal_coding') for variant in proximal_coding])
    sigma_m_list = list(map(sigma_m.get, df.most_severe_consequence.tolist()))
    df['sigma_m_var'] = sigma_m_list
    category_list = list(map(category_dict.get, df.most_severe_consequence.tolist()))
    df['category'] = category_list
    df = df.assign(sigma_m_1=1)
    df = df.assign(sigma_m_005=0.05)
    df = df[df.sigma_m_var.notnull()]
    return df

# Keep diagonals and multiples every other cell by .99
def is_pos_def(x):
    i = 0
    x = np.matrix(x)
    if np.all(np.linalg.eigvals(x) > 0):
        return x
    else:
        while not np.all(np.linalg.eigvals(x) > 0):
            x = 0.99*x + 0.01*np.diag(np.diag(x))
            break
            i += 1 
            if i >= 2:
                print("TAKING TOO LONG")
                break
        return x

def return_BF(U, beta, v_beta, mu):
    v_beta = is_pos_def(v_beta)
    v_beta_inv = np.linalg.inv(v_beta)
    U = is_pos_def(U)
    U_inv = np.linalg.inv(U)
    if v_beta_inv is not np.nan and U_inv is not np.nan:
        A2 = U_inv + v_beta_inv
        b2 = v_beta_inv * beta
        Abinv = np.linalg.lstsq(A2, b2, rcond=-1)[0]
        fat_middle = v_beta_inv - (v_beta_inv.dot(np.linalg.inv(U_inv + v_beta_inv))).dot(v_beta_inv)
        logBF = -0.5 * np.linalg.slogdet(npm.eye(beta.shape[0]) + v_beta_inv*U)[1] + 0.5 * beta.T.dot(v_beta_inv.dot(beta)) - 0.5*(((beta-mu).T).dot(fat_middle)).dot(beta-mu)
        logBF = float(np.array(logBF))
        log10BF = logBF/np.log(10)
        return log10BF

def assign_R_var(model, M):
    if model == 'independent_effects':
        R_var = np.diag(np.ones(M))
    elif model == 'similar_effects':
        R_var = np.ones((M, M))
    return R_var

def generate_beta_se_gene(disease_string, subset_df):
    se2_list = subset_df['SE'].tolist()
    if 'BETA' in subset_df.columns:
        beta_list = subset_df['BETA'].tolist()
    else: 
        beta_list = np.log(subset_df['OR'].tolist())
    return beta_list, se2_list

def calculate_all_params(df, disease_string, M, key, sigma_m_type, j, S):
    R_study = np.diag(np.ones(S))
    subset_df = df[df['gene_symbol'] == key]
    sigma_m = subset_df[sigma_m_type].tolist()
    diag_sigma_m = np.diag(np.atleast_1d(np.array(sigma_m)))
    R_var_indep = assign_R_var('independent_effects', M)
    R_var_sim = assign_R_var('similar_effects', M)
    S_var_indep = np.dot(np.dot(diag_sigma_m, R_var_indep), diag_sigma_m)
    S_var_sim = np.dot(np.dot(diag_sigma_m, R_var_sim), diag_sigma_m)
    beta_list, se2_list = generate_beta_se_gene(disease_string, subset_df)
    beta = np.array(beta_list).reshape(-1,1)
    mu = np.zeros(beta.shape)
    v_beta = np.diag(np.array(se2_list).reshape(-1))
    U_indep = np.kron(R_study, np.kron(S_var_indep, R_phen))
    U_sim = np.kron(R_study, np.kron(S_var_sim, R_phen))
    return U_indep, U_sim, beta, v_beta, mu

def run_mrp_gene_level(df, disease_string, S):
    m_dict = df.groupby('gene_symbol').size()
    bf_dict = {}
    bf_dfs = []
    sigma_m_types = ['sigma_m_var', 'sigma_m_005', 'sigma_m_1']
    for sigma_m_type in sigma_m_types:
        print('Sigma m type ' + sigma_m_type + ':')
        for i, (key, value) in enumerate(m_dict.items()):
            if i % 1000 == 0:
                print('Done ' + str(i) + ' genes out of ' + str(len(m_dict)))
            M = value
            U_indep, U_sim, beta, v_beta, mu = calculate_all_params(df, disease_string, M, key, sigma_m_type, i, S)
            bf_indep = return_BF(U_indep, beta, v_beta, mu)
            bf_sim = return_BF(U_sim, beta, v_beta, mu)
            bf_dict[key] = [bf_indep, bf_sim]
        bf_df = pd.DataFrame.from_dict(bf_dict, orient='index').reset_index().rename(columns={'index': 'gene_symbol', 0: 'log_10_BF_' + sigma_m_type + '_indep', 1: 'log_10_BF_' + sigma_m_type + '_sim'})
        bf_dfs.append(bf_df)
    inner_merge = partial(pd.merge, on='gene_symbol', how='inner')
    df = reduce(inner_merge, bf_dfs)
    return df

def run_mrp(df, disease_string, S):
    print('Running MRP for ' + disease_string + '...')
    df = run_mrp_gene_level(df, disease_string, S)
    return df

if __name__ == '__main__':
    disease_string = sys.argv[1]
    print("Reading in summary stats for " + disease_string + "...")
    #disease_df = read_in_summary_stats(disease_string)
    #disease_df = disease_df[disease_df['SE'].notnull()]
  
    #REMOVE THIS BLOCK IS TEST
    disease_df = pd.read_csv('head.tsv', sep='\t')
    disease_df.insert(loc=0, column='V', value=disease_df['#CHROM'].astype(str).str.cat(disease_df['POS'].astype(str), sep=':').str.cat(disease_df['REF'], sep=':').str.cat(disease_df['ALT'], sep=':'))
    disease_df = disease_df[disease_df['SE'].notnull()]

    print("Reading in metadata file...")
    metadata = read_metadata('/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-gene_consequence_wb_maf.tsv')
    merged = disease_df.merge(metadata)
    print("Before MAF filter:")
    print(len(merged))
    print("After MAF filter:")
    merged = merged[(merged.wb_maf <= 0.01) & (merged.wb_maf > 0)]
    print(len(merged))
    merged = set_sigmas(merged)
    S = 1
    K = 1 #looking at each phenotype separately, since otherwise, there will be correlated errors in our case. Can change if phenotypes are sufficiently different
    merged[['V', 'category']].to_csv('mrp_' + disease_string + '_categories.tsv', sep='\t', index=False)
    R_phen = np.diag(np.ones(K))
    result = run_mrp(merged, disease_string, S)
    result.sort_values('log_10_BF_sigma_m_var_indep').to_csv('mrp_' + disease_string + '_proximal_coding.tsv', sep='\t', index=False)
    for filter_out in ['proximal_coding', 'pav']:
        merged = merged[merged.category != filter_out]
        result = run_mrp(merged, disease_string, S)
        if filter_out == 'proximal_coding':
            result.sort_values('log_10_BF_sigma_m_var_indep').to_csv('mrp_' + disease_string + '_pav.tsv', sep='\t', index=False)
        else:
            result.sort_values('log_10_BF_sigma_m_var_indep').to_csv('mrp_' + disease_string + '_ptv.tsv', sep='\t', index=False)
