import os
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr, spearmanr
import sys

def set_up_df(phen_PATH, metadata, pheno):
    
    """
    Sets up dataframe for a phenotype for calculating Pearson
        correlations.

    Parameters:
    phen_PATH: Path to summary statistics for the phenotype.
    metadata: Path to metadata file containing LD independence information,
        MAF, etc.
    pheno: Phenotype name.

    Returns:
    df: Ready dataframe for Pearson correlation calculation.
    
    """
    if len(phen_PATH) != 0:
        df = pd.read_csv(
            phen_PATH,
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
    
    """
    Extends correlation and p-value lists with nan and 1 respectively if
        insufficient or bad data for pheno2.
    
    Parameters:
    p_corrs_nohla: Empirical estimates of Pearson correlation across phenotypes
        for the comparator pheno1.
    p_ps_nohla: Corresponding p-values.

    Returns:
    p_corrs_nohla: Empirical estimates of Pearson correlation across phenotypes
        for the comparator pheno1, + nan.
    p_ps_nohla: Corresponding p-values, + 1.

    """
    p_corrs_nohla.append(np.nan)
    p_ps_nohla.append(1)
    return p_corrs_nohla, p_ps_nohla

    
def build_R_phen_row(pheno1, phen_info, row_num):

    """
    Builds R_phen row for clustering (array data).

    Parameters:
    pheno1: Comparator phenotype.
    phen_info: phenotype_info.tsv file.
    row_num: Row in R_phen.

    Returns:
    p_corrs_nohla: Empirical estimates of Pearson correlation across phenotypes
        for the comparator pheno1.
    p_ps_nohla: Corresponding p-values.
    lens_nohla: Length of dataframe containing variants that meet the criteria:
        - SE <= 0.5
        - MAF >- 0.01
        - P (pheno1 or pheno2) <= 1e-5
        - LD independent

    """
    # Filter for SE as you read it in
    print("Reading in metadata...")
    metadata = pd.read_table('/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep_mpc_pli.tsv')
    phen1_PATH = list(phen_info[phen_info['GBE_ID'] == pheno1]['PATH'])[0]
    print("Filtering on p-value and SE...")
    df1 = set_up_df(phen1_PATH, metadata, pheno1)
    phenos = list(phen_info['GBE_ID'])
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
                phen2_PATH = list(phen_info[phen_info['GBE_ID'] == pheno2]['PATH'])[0]
                df2 = set_up_df(phen2_PATH, metadata, pheno2)
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
