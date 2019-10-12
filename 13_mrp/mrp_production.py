from __future__ import division
import pandas as pd
from functools import partial, reduce
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import numpy.matlib as npm
from numpy.linalg import LinAlgError
from scipy.stats import describe
from scipy.stats import beta
from collections import defaultdict
import subprocess
import math
import os
import sys
import argparse

# Read in summary statistics from GBE sumstats
def read_in_summary_stats(pops, phenos, datasets):
    sumstats_files = []
    for pop in pops:
        for pheno in phenos:
            for dataset in datasets:
                if dataset == "exome":
                    findCMD = (
                        "find /oak/stanford/groups/mrivas/ukbb24983/"
                        + dataset
                        + '/gwas/ -name "*'
                        + pheno
                        + '.*" | grep -v "exome-spb.log" | grep -v freeze | grep -v old | grep '
                        + pop
                    )
                elif dataset == "cal":
                    findCMD = (
                        "find /oak/stanford/groups/mrivas/ukbb24983/"
                        + dataset
                        + '/gwas/ -name "*'
                        + pheno
                        + '.*" | grep -v "genotyped.log" | grep -v freeze | grep -v old | grep '
                        + pop
                    )
                out = subprocess.Popen(
                    findCMD,
                    shell=True,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                (stdout, stderr) = out.communicate()
                sumstat_file = stdout.strip().decode("utf-8")
                print("Sumstat file:")
                print(sumstat_file)
                df = pd.read_csv(sumstat_file, sep="\t")
                df.insert(
                    loc=0,
                    column="V",
                    value=df["#CHROM"]
                    .astype(str)
                    .str.cat(df["POS"].astype(str), sep=":")
                    .str.cat(df["REF"], sep=":")
                    .str.cat(df["ALT"], sep=":"),
                )
                #sumstat_files.append(df)
                sumstat_files.append(sumstat_file)
    return sumstat_files

# Need to have: gene, consequence, pop/global allele frequencies.
def read_metadata(metadata_path):
    metadata = pd.read_csv(metadata_path, sep="\t")
    return metadata


def set_sigmas(df):
    ptv = [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "start_lost",
        "stop_lost",
    ]
    pav = [
        "protein_altering_variant",
        "inframe_deletion",
        "inframe_insertion",
        "splice_region_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "missense_variant",
    ]
    proximal_coding = [
        "synonymous_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "coding_sequence_variant",
        "incomplete_terminal_codon_variant",
        "TF_binding_site_variant",
    ]
    to_filter = [
        "regulatory_region_variant",
        "intron_variant",
        "intergenic_variant",
        "downstream_gene_variant",
        "mature_miRNA_variant",
        "non_coding_transcript_exon_variant",
        "upstream_gene_variant",
        "NA",
        "NMD_transcript_variant",
    ]
    print("Before consequence filter:")
    print(len(df))
    df = df[~df.most_severe_consequence.isin(to_filter)]
    print("After consequence filter:")
    print(len(df))
    print("Setting sigmas...")
    sigma_m_ptv = 0.2
    sigma_m_pav = 0.05
    sigma_m_pc = 0.05
    sigma_m = dict(
        [(variant, sigma_m_ptv) for variant in ptv]
        + [(variant, sigma_m_pav) for variant in pav]
        + [(variant, sigma_m_pc) for variant in proximal_coding]
    )
    category_dict = dict(
        [(variant, "ptv") for variant in ptv]
        + [(variant, "pav") for variant in pav]
        + [(variant, "proximal_coding") for variant in proximal_coding]
    )
    sigma_m_list = list(map(sigma_m.get, df.most_severe_consequence.tolist()))
    df["sigma_m_var"] = sigma_m_list
    category_list = list(map(category_dict.get, df.most_severe_consequence.tolist()))
    df["category"] = category_list
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
            x = 0.99 * x + 0.01 * np.diag(np.diag(x))
            break
            i += 1
            if i >= 2:
                print("TAKING TOO LONG")
                break
        return x


def safe_inv(X, matrix_name, gene):
    try:
        X_inv = np.linalg.inv(X)
    except LinAlgError as err:
        print("Could not invert " + matrix_name + " for gene " + gene + ":")
        print(X)
        return np.nan
    return X_inv


def return_BF(U, beta, v_beta, mu, gene):
    v_beta = is_pos_def(v_beta)
    v_beta_inv = safe_inv(v_beta, "v_beta", gene)
    U = is_pos_def(U)
    U_inv = safe_inv(U, "U", gene)
    if v_beta_inv is not np.nan and U_inv is not np.nan:
        A2 = U_inv + v_beta_inv
        b2 = v_beta_inv * beta
        try:
            Abinv = np.linalg.lstsq(A2, b2, rcond=-1)[0]
        except LinAlgError as err:
            return np.nan
        fat_middle = v_beta_inv - (
            v_beta_inv.dot(np.linalg.inv(U_inv + v_beta_inv))
        ).dot(v_beta_inv)
        logBF = (
            -0.5 * np.linalg.slogdet(npm.eye(beta.shape[0]) + v_beta_inv * U)[1]
            + 0.5 * beta.T.dot(v_beta_inv.dot(beta))
            - 0.5 * (((beta - mu).T).dot(fat_middle)).dot(beta - mu)
        )
        logBF = float(np.array(logBF))
        log10BF = logBF / np.log(10)
        return log10BF
    else:
        return np.nan


def assign_R_var(model, M):
    if model == "independent_effects":
        R_var = np.diag(np.ones(M))
    elif model == "similar_effects":
        R_var = np.ones((M, M))
    return R_var


def generate_beta_se_gene(subset_df):
    se2_list = subset_df["SE"].tolist()
    if "BETA" in subset_df.columns:
        beta_list = subset_df["BETA"].tolist()
    else:
        beta_list = np.log(subset_df["OR"].tolist())
    return beta_list, se2_list


def calculate_all_params(df, M, key, sigma_m_type, j, S):
    R_study = np.diag(np.ones(S))
    subset_df = df[df["gene_symbol"] == key]
    sigma_m = subset_df[sigma_m_type].tolist()
    diag_sigma_m = np.diag(np.atleast_1d(np.array(sigma_m)))
    R_var_indep = assign_R_var("independent_effects", M)
    R_var_sim = assign_R_var("similar_effects", M)
    S_var_indep = np.dot(np.dot(diag_sigma_m, R_var_indep), diag_sigma_m)
    S_var_sim = np.dot(np.dot(diag_sigma_m, R_var_sim), diag_sigma_m)
    beta_list, se2_list = generate_beta_se_gene(subset_df)
    beta = np.array(beta_list).reshape(-1, 1)
    mu = np.zeros(beta.shape)
    v_beta = np.diag(np.square(np.array(se2_list)).reshape(-1))
    U_indep = np.kron(R_study, np.kron(S_var_indep, R_phen))
    U_sim = np.kron(R_study, np.kron(S_var_sim, R_phen))
    return U_indep, U_sim, beta, v_beta, mu


def run_mrp_gene_level(df, S, mode):
    m_dict = df.groupby("gene_symbol").size()
    bf_dict = {}
    bf_dfs = []
    sigma_m_types = ["sigma_m_var", "sigma_m_005", "sigma_m_1"]
    for sigma_m_type in sigma_m_types:
        print("Sigma m type " + sigma_m_type + ":")
        for i, (key, value) in enumerate(m_dict.items()):
            if i % 1000 == 0:
                print("Done " + str(i) + " genes out of " + str(len(m_dict)))
            M = value
            U_indep, U_sim, beta, v_beta, mu = calculate_all_params(
                df, M, key, sigma_m_type, i, S
            )
            bf_indep = return_BF(U_indep, beta, v_beta, mu, key)
            bf_sim = return_BF(U_sim, beta, v_beta, mu, key)
            bf_dict[key] = [bf_indep, bf_sim]
        bf_df = (
            pd.DataFrame.from_dict(bf_dict, orient="index")
            .reset_index()
            .rename(
                columns={
                    "index": "gene_symbol",
                    0: "log_10_BF_" + sigma_m_type + "_indep_" + mode,
                    1: "log_10_BF_" + sigma_m_type + "_sim_" + mode,
                }
            )
        )
        bf_dfs.append(bf_df)
    inner_merge = partial(pd.merge, on="gene_symbol", how="inner")
    df = reduce(inner_merge, bf_dfs)
    return df


def run_mrp(df, disease_string, S, mode):
    df = run_mrp_gene_level(df, disease_string, S, mode)
    return df


def merge_and_filter(pops, phenos, variant_filter, datasets):

    print("Reading in summary stats for:")
    print("Populations: " + ", ".join(pops))
    print("Phenotypes: "  + ", ".join(phenos))
    print("Datasets: " + ", ".join(datasets))

    disease_df, sumstat_file = read_in_summary_stats(pops, phenos, datasets)
    print("Before SE filter...")
    print(len(disease_df))
    disease_df = disease_df[disease_df["SE"].notnull()]
    disease_df = disease_df[disease_df["SE"] <= 0.5]
    print("After SE filter...")
    print(len(disease_df))
    mrp_prefix = os.path.dirname(sumstat_file.replace("gwas", "mrp"))

    # disease_df = pd.read_csv('head.tsv', sep='\t')
    # disease_df = disease_df[disease_df['SE'].notnull()]
    # disease_df.insert(loc=0, column='V', value=disease_df['#CHROM'].astype(str).str.cat(disease_df['POS'].astype(str), sep=':').str.cat(disease_df['REF'], sep=':').str.cat(disease_df['ALT'], sep=':'))
    # mrp_prefix=""

    # Merge metadata
    print("Reading in metadata file...")
    if mode == "exome":
        metadata = read_metadata(
            "/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-gene_consequence_wb_maf_final.tsv"
        )
    elif mode == "cal":
        metadata = read_metadata(
            "/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-gene_consequence_wb_maf_final.tsv"
        )
    merged = disease_df.merge(metadata)
    print("Before MAF filter:")
    print(len(merged))
    print("After MAF filter:")
    merged = merged[(merged.wb_maf <= 0.01) & (merged.wb_maf > 0)]
    print(len(merged))
    merged = set_sigmas(merged)
    return merged, mrp_prefix


def generate_results(merged, disease_string, S):
    results = []
    print("Running proximal coding variants...")
    results.append(run_mrp(merged, disease_string, S, "proximal_coding"))
    for filter_out in ["proximal_coding", "pav"]:
        merged = merged[merged.category != filter_out]
        if filter_out == "proximal_coding":
            print("Running PAVs...")
            results.append(run_mrp(merged, disease_string, S, "pav"))
        else:
            print("Running PTVs...")
            results.append(run_mrp(merged, disease_string, S, "ptv"))
    inner_merge = partial(pd.merge, on="gene_symbol", how="outer")
    df = reduce(inner_merge, results)
    return df


def return_input_args(args):
    S = len(args.pops)
    K = len(args.phenos)
    R_study = [
        np.diag(np.ones(S)) if x == "independent" else np.ones((S, S))
        for x in args.R_study_model
    ]
    R_phen = [
        np.diag(np.ones(K)) if x == "independent" else np.ones((K, K))
        for x in args.R_phen_model
    ]
    return (
        S,
        K,
        args.pops,
        args.phenos,
        R_study,
        R_phen,
        args.agg,
        args.sigma_type,
        args.variant_filter,
        args.datasets,
    )


if __name__ == "__main__":
    # Argparse should get all of the parameters needed:
    # S, M, and K should be specified
    parser = argparse.ArgumentParser(
        description="MRP takes in several variables that affect how it runs."
    )
    parser.add_argument(
        "--S",
        choices=["white_british", "non_british_white", "african", "e_asian", "s_asian"],
        type=str,
        nargs="+",
        required=True,
        dest="pops",
        help="names of populations/studies to be meta-analyzed; at least one required",
    )
    parser.add_argument(
        "--R_study",
        choices=["independent", "similar"],
        type=str,
        nargs="+",
        default=["independent"],
        dest="R_study_model",
        help="type of model across studies. options: independent, similar (default: independent). can run both.",
    )
    # ERROR CHECKING FOR ELIGIBLE PHENOS?
    parser.add_argument(
        "--K",
        type=str,
        nargs="+",
        required=True,
        dest="phenos",
        help="names of phenotypes to be studied jointly; at least one required",
    )
    parser.add_argument(
        "--R_phen",
        choices=["independent", "similar"],
        type=str,
        nargs="+",
        default=["independent"],
        dest="R_phen_model",
        help="type of model across phenotypes. options: independent, similar (default: independent). can run both.",
    )
    parser.add_argument(
        "--M",
        choices=["variant", "gene"],
        type=str,
        nargs="+",
        default=["gene"],
        dest="agg",
        help="unit of aggregation. options: variant, gene (default: gene). can run both.",
    )
    parser.add_argument(
        "--sigma_type",
        choices=["var", "1", "0.05"],
        type=str,
        nargs="+",
        default=["var"],
        dest="sigma_type",
        help="scaling factor for variants. options: var (i.e. 0.2 for ptvs, 0.05 for pavs/pcvs), 1, 0.05 (default: var). can run multiple.",
    )
    parser.add_argument(
        "--variants",
        choices=["pcv", "pav", "ptv"],
        type=str,
        nargs="+",
        default=["ptv"],
        dest="variant_filter",
        help="variants to consider. options: proximal coding [pcv], protein-altering [pav], protein truncating [ptv] (default: ptv). can run multiple.",
    )
    parser.add_argument(
        "--datasets",
        choices=["cal", "exome"],
        type=str,
        nargs="+",
        default=["cal"],
        dest="datasets",
        help="which UKBB dataset to use. options: cal, exome (default: cal). can run multiple.",
    )
    args = parser.parse_args()

    S, K, pops, phenos, R_study, R_phen, agg, sigma_type, variant_filter, datasets = return_input_args(
        args
    )
    
    merged, mrp_prefix = merge_and_filter(pops, phenos, variant_filter, datasets)
    # S = 1
    # K = 1 #looking at each phenotype separately, since otherwise, there will be correlated errors in our case. Can change if phenotypes are sufficiently different
    # gene_pos = merged[['gene_symbol', '#CHROM', 'POS']]
    # idx = gene_pos.groupby(['gene_symbol'])['POS'].transform(min) == gene_pos['POS']
    # gene_pos = gene_pos[idx]
    # merged[['V', 'gene_symbol', 'category']].to_csv(os.path.join(mrp_prefix, disease_string + '_categories.tsv'), sep='\t', index=False)
    # R_phen = np.diag(np.ones(K))
    # print('Running MRP for ' + disease_string + '...')
    # df = generate_results(merged, disease_string, S)
    # df = df.merge(gene_pos, on='gene_symbol', how='inner')
    # df.to_csv(os.path.join(mrp_prefix, disease_string + '_gene.tsv'), sep='\t', index=False)
