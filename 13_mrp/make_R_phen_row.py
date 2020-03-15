import os
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr, spearmanr
import sys

def set_up_df(phen_path, metadata, pheno):
    if len(phen_path) != 0:
        df = pd.read_csv(
            phen_path,
            sep="\t",
            dtype={
                "#CHROM": str,
                "POS": np.int32,
                "ID": str,
                "REF": str,
                "ALT": str,
                "A1": str,
                "FIRTH?": str,
                "TEST": str,
            },
        )
        df.insert(
            loc=0,
            column="V",
            value=df["#CHROM"]
            .astype(str)
            .str.cat(df["POS"].astype(str), sep=":")
            .str.cat(df["REF"], sep=":")
            .str.cat(df["ALT"], sep=":"),
        )
        df = df[df["SE"].notnull()]
        df = df[df["SE"].astype(float) <= 0.5]
        if "OR" in df.columns:
            df["BETA"] = np.log(df["OR"].astype("float64"))
        df = df.merge(metadata)
        df = df[(df.maf >= 0.01) & (df.ld_indep == True)]
        df = df[["V", "BETA", "P"]]
        df.columns = ["V", "BETA_" + pheno, "P_" + pheno]
        return df
    else:
        return []

def extend_lists_w_missing(p_corrs_nohla, p_ps_nohla):
    p_corrs_nohla.append(np.nan)
    p_ps_nohla.append(1)
    return p_corrs_nohla, p_ps_nohla

    
def build_R_phen_row(pheno1, phen_info, row_num):

    """
    Builds R_phen using phen_corr (calculated using the method directly above this).

    Parameters:
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.
    map_file: Input file containing summary statistic paths + pop and pheno data.

    Returns:
    R_phen: Empirical estimates of genetic correlation across phenotypes.

    """
    # Filter for SE as you read it in
    print("Reading in metadata...")
    metadata = pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv')
    phen1_path = list(phen_info[phen_info['PHEN'] == pheno1]['path'])[0]
    print("Filtering on p-value and SE...")
    df1 = set_up_df(phen1_path, metadata, pheno1)
    phenos = list(phen_info['PHEN'])
    if len(df1) == 0:
        p_corrs_nohla = np.repeat(np.nan, len(phenos), axis=0)
        p_ps_nohla = np.repeat(np.nan, len(phenos), axis=0)
        p_corrs_nohla[int(row_num) - 1] = 1
        p_ps_nohla[int(row_num) - 1] = 0
        lens_nohla = np.repeat(0, len(phenos), axis=0)
        return p_corrs_nohla, p_ps_nohla, lens_nohla
    p_corrs_nohla, p_ps_nohla, lens_nohla = [], [], []
    compute = False
    for k2, pheno2 in enumerate(phenos):
        if pheno1 == pheno2:
            p_corrs_nohla.append(1)
            p_ps_nohla.append(0)
            df = df1[df1["P_" + pheno1] <= 1e-5]
            # FILTER OUT HLA
            df = df[~((df["V"].str.split(":").apply(lambda x: x[0]).astype(int) == 6) & (df["V"].str.split(":").apply(lambda x: x[1]).astype(int).between(25477797, 36448354)))]
            lens_nohla.append(len(df))
            compute = True
        else:
            if compute == True:
                phen2_path = list(phen_info[phen_info['PHEN'] == pheno2]['path'])[0]
                df2 = set_up_df(phen2_path, metadata, pheno2)
                if len(df2) != 0:
                    # FILTER OUT HLA
                    df = df1.merge(df2)
                    df = df[~((df["V"].str.split(":").apply(lambda x: x[0]).astype(int) == 6) & (df["V"].str.split(":").apply(lambda x: x[1]).astype(int).between(25477797, 36448354)))]
                    df = df[(df["P_" + pheno1] <= 1e-5) | (df["P_" + pheno2] <= 1e-5)]
                    beta1, beta2 = list(df["BETA_" + pheno1]), list(df["BETA_" + pheno2])
                    lens_nohla.append(len(beta1))
                    if ((len(beta1) > 2) and (len(beta2) > 2)):
                        p_corr_nohla, p_p_nohla = pearsonr(beta1, beta2)
                        p_corrs_nohla.append(p_corr_nohla)
                        p_ps_nohla.append(p_p_nohla)
                    else:
                        p_corrs_nohla, p_ps_nohla = extend_lists_w_missing(p_corrs_nohla, p_ps_nohla)
                else:
                    p_corrs_nohla, p_ps_nohla = extend_lists_w_missing(p_corrs_nohla, p_ps_nohla)
                    lens_nohla.append(0)
            else:
                p_corrs_nohla, p_ps_nohla = extend_lists_w_missing(p_corrs_nohla, p_ps_nohla)
                lens_nohla.append(0)
    return p_corrs_nohla, p_ps_nohla, lens_nohla


if __name__ == "__main__":
    pheno1 = sys.argv[1]
    row_num = sys.argv[2]
    phen_info = pd.read_table("sumstat_paths.tsv")
    phen_info = phen_info.fillna("")
    p_corrs_nohla, p_ps_nohla, lens_nohla = build_R_phen_row(pheno1, phen_info, row_num)
    if (len(p_corrs_nohla) == len(p_ps_nohla)) and (len(p_ps_nohla) == len(lens_nohla)):
        np.save("R_phen_rows/" + str(row_num) + "_pearson_nohla.npy", p_corrs_nohla)
        np.save("R_phen_rows/" + str(row_num) + "_pearson_p_nohla.npy", p_ps_nohla) 
        np.save("R_phen_rows/" + str(row_num) + "_betalen_nohla.npy", lens_nohla)
    else:
        raise ValueError("Lengths of all vectors not equal!")
