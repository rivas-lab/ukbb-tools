from __future__ import division
import argparse

# Read in summary statistics from GBE sumstats
def read_in_summary_stats(pops, phenos, datasets, conserved_columns):
    file_paths = []
    sumstat_files = []
    for pop in pops:
        for pheno in phenos:
            for dataset in datasets:
                findCMD = (
                    "find /oak/stanford/groups/mrivas/ukbb24983/"
                    + dataset
                    + '/gwas/ -name "*.gz" | grep '
                    + pop
                    + " | grep -v freeze | grep -v old | grep "
                    + pheno
                )
                out = subprocess.Popen(
                    findCMD,
                    shell=True,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                (stdout, stderr) = out.communicate()
                file_path = stdout.strip().decode("utf-8")
                if os.path.exists(file_path):
                    print(file_path)
                    df = pd.read_csv(file_path, sep="\t")
                    df.insert(
                        loc=0,
                        column="V",
                        value=df["#CHROM"]
                        .astype(str)
                        .str.cat(df["POS"].astype(str), sep=":")
                        .str.cat(df["REF"], sep=":")
                        .str.cat(df["ALT"], sep=":"),
                    )
                    # Filter for SE as you read it in
                    df = rename_columns(df, conserved_columns, pop, pheno)
                    df = df[df["SE" + "_" + pop + "_" + pheno].notnull()]
                    df = df[df["SE" + "_" + pop + "_" + pheno] <= 0.5]
                    sumstat_files.append(df)
                    file_paths.append(file_path)
                else:
                    print(
                        "A summary statistic file cannot be found for population: {}; phenotype: {}; dataset: {}. A placeholder df will be added instead.".format(
                            pop, pheno, dataset
                        )
                    )
                    if ("QT" in pheno) or ("INI" in pheno):
                        df = pd.DataFrame(
                            columns=conserved_columns
                            + ["TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P"]
                        )
                    else:
                        df = pd.DataFrame(
                            columns=conserved_columns
                            + ["FIRTH?", "TEST", "OBS_CT", "OR", "SE", "Z_STAT", "P"]
                        )
                    df = rename_columns(df, conserved_columns, pop, pheno)
                    sumstat_files.append(df)
                    file_paths.append("file_not_found")
    return file_paths, sumstat_files


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

    df = df[~df.most_severe_consequence.isin(to_filter)]
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


def delete_rows_and_columns(matrix, indices_to_remove):
    matrix = np.delete(matrix, indices_to_remove, axis=0)
    matrix = np.delete(matrix, indices_to_remove, axis=1)
    return matrix


def adjust_U_for_missingness(U, v_beta, beta, beta_list):
    indices_to_remove = np.argwhere(np.isnan(beta_list))
    U = delete_rows_and_columns(U, indices_to_remove)
    v_beta = delete_rows_and_columns(v_beta, indices_to_remove)
    beta = beta[~np.isnan(beta)].reshape(-1, 1)
    return U, v_beta, beta


def generate_beta_se(subset_df, pops, phenos):
    beta_list = []
    se_list = []
    for pop in pops:
        for index, row in subset_df.iterrows():
            for pheno in phenos:
                if "BETA" + "_" + pop + "_" + pheno in subset_df.columns:
                    beta_list.append(row["BETA" + "_" + pop + "_" + pheno])
                elif "OR" + "_" + pop + "_" + pheno in subset_df.columns:
                    beta_list.append(np.log(row["OR" + "_" + pop + "_" + pheno]))
                se_list.append(row["SE" + "_" + pop + "_" + pheno])
    return beta_list, se_list


def calculate_all_params(
    df, pops, phenos, M, key, sigma_m_type, j, S, R_study, R_phen, R_var_model
):
    subset_df = df[df["gene_symbol"] == key]
    sigma_m = subset_df[sigma_m_type].tolist()
    diag_sigma_m = np.diag(np.atleast_1d(np.array(sigma_m)))
    R_var = np.diag(np.ones(M)) if R_var_model == "independent" else np.ones((M, M))
    S_var = np.dot(np.dot(diag_sigma_m, R_var), diag_sigma_m)
    beta_list, se_list = generate_beta_se(subset_df, pops, phenos)
    beta = np.array(beta_list).reshape(-1, 1)
    v_beta = np.diag(np.square(np.array(se_list)).reshape(-1))
    U = np.kron(R_study, np.kron(S_var, R_phen))
    U, v_beta, beta = adjust_U_for_missingness(U, v_beta, beta, beta_list)
    mu = np.zeros(beta.shape)
    return U, beta, v_beta, mu


def run_mrp(
    dfs,
    S,
    K,
    pops,
    phenos,
    R_study,
    R_study_model,
    R_phen,
    R_phen_model,
    R_var_model,
    analysis,
    sigma_m_type,
    conserved_columns,
):
    outer_merge = partial(pd.merge, on=conserved_columns, how="outer")
    df = reduce(outer_merge, dfs)
    m_dict = df.groupby("gene_symbol").size()
    bf_dict = {}
    for i, (key, value) in enumerate(m_dict.items()):
        if i % 1000 == 0:
            print("Done " + str(i) + " genes out of " + str(len(m_dict)))
        M = value
        U, beta, v_beta, mu = calculate_all_params(
            df, pops, phenos, M, key, sigma_m_type, i, S, R_study, R_phen, R_var_model
        )
        bf = return_BF(U, beta, v_beta, mu, key)
        bf_dict[key] = bf
    bf_df = (
        pd.DataFrame.from_dict(bf_dict, orient="index")
        .reset_index()
        .rename(
            columns={
                "index": "gene_symbol",
                0: "log_10_BF"
                + "_study_"
                + R_study_model
                + "_phen_"
                + R_phen_model
                + "_var_"
                + R_var_model
                + "_"
                + sigma_m_type
                + "_"
                + analysis,
            }
        )
    )
    return bf_df


def filter_category(sumstats_files, variant_filter):
    analysis_files = []
    for df in sumstats_files:
        if variant_filter == "ptv":
            df = df[df.category == "ptv"]
        elif variant_filter == "pav":
            df = df[(df.category == "ptv") | (df.category == "pav")]
        analysis_files.append(df)
    return analysis_files


def rename_columns(df, conserved_columns, pop, pheno):
    if "LOG(OR)_SE" in df.columns:
        df.rename(columns={"LOG(OR)_SE": "SE"}, inplace=True)
    columns_to_rename = list(set(df.columns) - set(conserved_columns))
    renamed_columns = [(x + "_" + pop + "_" + pheno) for x in columns_to_rename]
    df.rename(columns=dict(zip(columns_to_rename, renamed_columns)), inplace=True)
    return df


def collect_and_filter(pops, phenos, datasets, conserved_columns):
    print("Reading in summary stats for:")
    print("Populations: " + ", ".join(pops))
    print("Phenotypes: " + ", ".join(phenos))
    print("Datasets: " + ", ".join(datasets))
    print("")

    file_paths, sumstat_files = read_in_summary_stats(
        pops, phenos, datasets, conserved_columns
    )

    filtered_sumstat_files = []
    exome_metadata = read_metadata(
        "/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-gene_consequence_wb_maf_final.tsv"
    )
    cal_metadata = read_metadata(
        "/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-gene_consequence_wb_maf_final.tsv"
    )
    print("")
    print("Filtering sumstats on MAF, setting sigmas, and filtering on consequence...")
    print("")
    for file_path, sumstat_file in zip(file_paths, sumstat_files):
        # Merge metadata
        if "exome" in file_path:
            sumstat_file = sumstat_file.merge(exome_metadata)
        elif "cal" in file_path:
            sumstat_file = sumstat_file.merge(cal_metadata)
        sumstat_file = sumstat_file[
            (sumstat_file.wb_maf <= 0.01) & (sumstat_file.wb_maf > 0)
        ]
        sumstat_file = set_sigmas(sumstat_file)
        filtered_sumstat_files.append(sumstat_file)
    return filtered_sumstat_files


def return_input_args(args):
    for arg in vars(args):
        setattr(args, arg, list(set(getattr(args, arg))))
    S = len(args.pops)
    K = len(args.phenos)
    R_study = [
        np.diag(np.ones(S)) if x == "independent" else np.ones((S, S))
        for x in args.R_study_models
    ]
    R_phen = [
        np.diag(np.ones(K)) if x == "independent" else np.ones((K, K))
        for x in args.R_phen_models
    ]
    return (
        S,
        K,
        args.pops,
        args.phenos,
        R_study,
        args.R_study_models,
        R_phen,
        args.R_phen_models,
        args.R_var_models,
        args.agg,
        args.sigma_m_types,
        args.variants,
        args.datasets,
    )


def print_params(
    analysis, R_study_model, R_phen_model, R_var_model, agg_type, sigma_m_type
):
    print("Analysis: " + analysis)
    print("R_study model: " + R_study_model)
    print("R_phen model: " + R_phen_model)
    print("R_var model: " + R_var_model)
    print("Aggregation by: " + agg_type)
    print("Variant weighting factor: " + sigma_m_type)
    print("")

if __name__ == "__main__":
    with open("../05_gbe/phenotype_info.tsv", "r") as phe_file:
        valid_phenos = [line.split()[0] for line in phe_file][1:]

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
        dest="R_study_models",
        help="type of model across studies. options: independent, similar (default: independent). can run both.",
    )
    parser.add_argument(
        "--K",
        choices=valid_phenos,
        metavar="PHENO1 PHENO2 ...",
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
        dest="R_phen_models",
        help="type of model across phenotypes. options: independent, similar (default: independent). can run both.",
    )
    parser.add_argument(
        "--R_var",
        choices=["independent", "similar"],
        type=str,
        nargs="+",
        default=["independent"],
        dest="R_var_models",
        help="type of model across variants. options: independent, similar (default: independent). can run both.",
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
        "--sigma_m_types",
        choices=["sigma_m_var", "sigma_m_1", "sigma_m_0.05"],
        type=str,
        nargs="+",
        default=["sigma_m_var"],
        dest="sigma_m_types",
        help="scaling factor for variants. options: var (i.e. 0.2 for ptvs, 0.05 for pavs/pcvs), 1, 0.05 (default: var). can run multiple.",
    )
    parser.add_argument(
        "--variants",
        choices=["pcv", "pav", "ptv"],
        type=str,
        nargs="+",
        default=["ptv"],
        dest="variants",
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
    print("Valid arguments. Importing required packages...")
    print("")
    import pandas as pd
    from functools import partial, reduce

    pd.options.mode.chained_assignment = None
    import numpy as np
    import numpy.matlib as npm
    from numpy.linalg import LinAlgError
    from scipy.stats import describe, beta
    from collections import defaultdict
    import subprocess
    import math
    import os

    S, K, pops, phenos, R_study_list, R_study_models, R_phen_list, R_phen_models, R_var_models, agg, sigma_m_types, variant_filters, datasets = return_input_args(
        args
    )

    conserved_columns = [
        "V",
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "A1",
        "most_severe_consequence",
        "wb_maf",
        "gene_symbol",
        "sigma_m_var",
        "category",
        "sigma_m_1",
        "sigma_m_005",
    ]
    sumstat_files = collect_and_filter(pops, phenos, datasets, conserved_columns)

    for dataset in datasets:
        for analysis in variant_filters:
            analysis_files = filter_category(sumstat_files, analysis)
            for R_study, R_study_model in zip(R_study_list, R_study_models):
                for R_phen, R_phen_model in zip(R_phen_list, R_phen_models):
                    for agg_type in agg:
                        for R_var_model in R_var_models:
                            for sigma_m_type in sigma_m_types:
                                print_params(
                                    analysis,
                                    R_study_model,
                                    R_phen_model,
                                    R_var_model,
                                    agg_type,
                                    sigma_m_type,
                                )
                                bf_df = run_mrp(
                                    analysis_files,
                                    S,
                                    K,
                                    pops,
                                    phenos,
                                    R_study,
                                    R_study_model,
                                    R_phen,
                                    R_phen_model,
                                    R_var_model,
                                    analysis,
                                    sigma_m_type,
                                    conserved_columns,
                                )
                                out_file = (
                                    "/oak/stanford/groups/mrivas/ukbb24983/"
                                    + dataset
                                    + "/mrp/"
                                    + "_".join(pops)
                                    + "_"
                                    + "_".join(phenos)
                                    + "_study_"
                                    + R_study_model
                                    + "_phen_"
                                    + R_phen_model
                                    + "_var_"
                                    + R_var_model
                                    + "_"
                                    + sigma_m_type
                                    + "_"
                                    + analysis
                                    + ".tsv"
                                )
                                bf_df.sort_values(
                                    by=list(set(bf_df.columns) - set(["gene_symbol"])),
                                    ascending=False,
                                ).to_csv(out_file, sep="\t", index=False)
                                print("Results written to " + out_file + ".")
