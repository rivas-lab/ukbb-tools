import pandas as pd
import numpy as np

def calculate_phen(a, b, pop1, pheno1, pop2, pheno2, df, pop_pheno_tuples):

    """
    Calculates a single entry in the phen_corr matrix.
    
    Parameters:
    pop1: Name of first population.
    pheno1: Name of first phenotype.
    pop2: Name of second population.
    pheno2: Name of second phenotype.
    df: Dataframe containing significant, common, LD-independent variants.
    pop_pheno_tuples: Indicate which populations/phenotypes to use to build R_phen.

    Returns:
    phen_corr[a, b]: One entry in the phen_corr matrix.
    """
    # If in lower triangle, do not compute; symmetric matrix
    if (a > b) or (a == b):
        return np.nan
    else:
        # if this combination of pop, pheno doesn't exist in the map file, then nan
        if ((pop1, pheno1) in pop_pheno_tuples) and (
            (pop2, pheno2) in pop_pheno_tuples
        ):
            
    beta1 = list(df["BETA_" + pop1 + "_" + pheno1])
    beta2 = list(df["BETA_" + pop2 + "_" + pheno2])

            phen_beta1, phen_beta2 = get_betas(df, pop1, pheno1, pop2, pheno2, "sig")
            return (
                pearsonr(phen_beta1, phen_beta2)[0]
                if phen_beta1 is not None
                else np.nan
            )
        else:
            return np.nan

def filter_for_phen_corr(df, map_file):

    """
    Filters the initial dataframe for the criteria used to build R_phen.

    Parameters:
    df: Merged dataframe containing all summary statistics.
    map_file: Dataframe indicating which summary statistics to use to build R_phen.

    Returns:
    df: Filtered dataframe that contains significant, common, LD-independent variants.

    """

    cols_to_keep = ["V", "maf", "ld_indep"]
    for col_type in "BETA_", "P_":
        cols_to_keep.extend(
            [col_type + pop + "_" + pheno for pop, pheno in pop_pheno_tuples]
        )
    df = df[cols_to_keep]
    # Get only LD-independent, common variants
    df = df[(df.maf >= 0.01) & (df.ld_indep == True)]
    df = df[
        (df["P_" + pop1 + "_" + pheno1] < 1e-5)
        | (df["P_" + pop2 + "_" + pheno2] < 1e-5)
    ]
    return df


def build_R_phen(phenos):

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
    K = len(phenos)
    R_phen = np.zeros((K, K))
    for k1, pheno1 in zip(range(K), phenos):
        for k2, pheno2 in zip(range(K), phenos):
            if k1 == k2:
                R_phen[k1, k2] = 1
            elif k1 > k2:
                R_phen[k1, k2] = R_phen[k2, k1]
            else:
                # Find sumstat files, load them in
                # Filter
                # Inner merge
                # get betas, calculate corrs
    R_phen = np.nan_to_num(R_phen)
    R_phen[abs(R_phen) < 0.01] = 0
    return R_phen


if __name__ == "__main__":
    phen_info = pd.read_table("../05_gbe/phenotype_info.tsv")
    phenos = list(phen_info[(phen_info["N"] >= 10) & ~(phen_info["#GBE_ID"].str.startswith("MED"))]["#GBE_ID"])
    R_phen = build_R_phen(phenos)
