from __future__ import division
import argparse


def read_in_summary_stats(pops, phenos, datasets, conserved_columns):

    """ 
    Reads in GBE summary statistics from the Rivas Lab file organization system on Sherlock.
  
    Additionally: adds a variant identifier ("V"), renames columns, and filters on SE (<= 0.5).
    Contains logic for handling the case that a summary statistic file is not found.
  
    Parameters: 
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    datasets: Unique list of datasets ("cal"/"exome") to use for analysis.
    conserved_columns: Columns in every annotated summary statistic file; basis of merges.
  
    Returns: 
    file_paths: List of strings containing file paths.
    sumstat_files: List of dataframes containing summary statistics.
  
    """

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


def set_sigmas(df):

    """ 
    Assigns appropriate sigmas to appropriate variants by annotation.
  
    Filters out variants not of interest;
    Sets sigmas by functional annotation;
    Additionally adds two extra columns for two other standard choices of a uniform sigma (1 and 0.05).
  
    Parameters: 
    df: Merged dataframe containing all variants across all studies and phenotypes.
  
    Returns: 
    df: Merged dataframe with three additional columns:
        sigma_m_var: Column of sigma values (mapped to functional annotation via the lists inside this method).
            NOTE: One can change the sigmas associated with each type of variant by adjusting this method.
        sigma_m_1: Uniform column of 1.
        sigma_m_005: Uniform column of 0.05.
  
    """

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


def is_pos_def(x):

    """ 
    Ensures a matrix is positive definite.
  
    Keep diagonals and multiples every other cell by .99
  
    Parameters: 
    x: Matrix to verify.
  
    Returns: 
    x: Verified (and, if applicable, adjusted) matrix.
  
    """

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


def safe_inv(X, matrix_name, block, agg_type):

    """ 
    Safely inverts a matrix, or returns NaN.
  
    Parameters: 
    X: Matrix to invert.
    matrix_name: One of "U"/"v_beta" - used to print messages when inversion fails.
    block: Name of the aggregation block (gene/variant) - used to print messages when inversion fails.
    agg_type: One of "gene"/"variant". Dictates block of aggregation. Used to print messages when inversion fails.
  
    Returns: 
    X_inv: Inverse of X.
  
    """

    try:
        X_inv = np.linalg.inv(X)
    except LinAlgError as err:
        print(
            "Could not invert " + matrix_name + " for " + agg_type + " " + block + "."
        )
        return np.nan
    return X_inv


def return_BF(U, beta, v_beta, mu, block, agg_type):

    """ 
    Given quantities calculated previously and the inputs, returns the associated Bayes Factor.
  
    Parameters: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K) dictating correlation structures; no missing data.
    beta: Effect size vector without missing data.
    v_beta: Diagonal matrix of variances of effect sizes without missing data.
    mu: A mean of genetic effects, size of beta (NOTE: default is 0, can change in the code below).
    block: Name of the aggregation block (gene/variant).
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
  
    Returns: 
    log10BF: log_10 Bayes Factor (ratio of marginal likelihoods of alternative model, which accounts for priors, and null).
  
    """

    v_beta = is_pos_def(v_beta)
    v_beta_inv = safe_inv(v_beta, "v_beta", block, agg_type)
    U = is_pos_def(U)
    U_inv = safe_inv(U, "U", block, agg_type)
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

    """ 
    Helper function to delete rows and columns from a matrix.
  
    Parameters: 
    matrix: Matrix that needs adjustment.
    indices_to_remove: Rows and columns to be deleted
  
    Returns: 
    matrix: Smaller matrix that has no missing data.
  
    """

    matrix = np.delete(matrix, indices_to_remove, axis=0)
    matrix = np.delete(matrix, indices_to_remove, axis=1)
    return matrix


def adjust_for_missingness(U, omega, beta, se, beta_list):
    """ 
    Deletes rows and columns where we do not have effect sizes/standard errors.

    Calls method delete_rows_and_columns, a helper function that calls the numpy command.
  
    Parameters: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K) dictating correlation structures; may relate to missing data.
    omega: (S*M*K x S*M*K) matrix that contains correlation of errors across variants, studies, and phenotypes. 
    beta: Vector of effect sizes within the unit of aggregation; may contain missing data.
    se: Vector of standard errors within the unit of aggregation; may contain missing data.
    beta_list: List of effect sizes within the unit of aggregation; may contain missing data.
  
    Returns: 
    U: Potentially smaller U matrix not associated with missing data.
    omega: Potentially smaller omega matrix not associated with missing data.
    beta: Potentially smaller beta vector without missing data.
    se: Potentially smaller SE vector without missing data.
  
    """
    indices_to_remove = np.argwhere(np.isnan(beta_list))
    U = delete_rows_and_columns(U, indices_to_remove)
    omega = delete_rows_and_columns(omega, indices_to_remove)
    beta = beta[~np.isnan(beta)].reshape(-1, 1)
    se = se[~np.isnan(se)]
    return U, omega, beta, se


def generate_beta_se(subset_df, pops, phenos):

    """ 
    Gathers effect sizes and standard errors from a unit of aggregation (gene/variant).
  
    Parameters: 
    subset_df: slice of the original df that encompasses the current unit of aggregation (gene/variant).
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
  
    Returns: 
    beta_list: A list of effect sizes (some may be missing) from the subset.
    se_list: A list of standard errors (some may be missing) from the subset.
  
    """

    beta_list = []
    se_list = []
    for pop in pops:
        for pheno in phenos:
            for index, row in subset_df.iterrows():
                if "BETA" + "_" + pop + "_" + pheno in subset_df.columns:
                    beta_list.append(row["BETA" + "_" + pop + "_" + pheno])
                elif "OR" + "_" + pop + "_" + pheno in subset_df.columns:
                    beta_list.append(np.log(row["OR" + "_" + pop + "_" + pheno]))
                se_list.append(row["SE" + "_" + pop + "_" + pheno])
    return beta_list, se_list


def calculate_all_params(
    df,
    pops,
    phenos,
    key,
    sigma_m_type,
    R_study,
    R_phen,
    R_var_model,
    agg_type,
    M,
    err_corr,
):

    """ 
    Calculates quantities needed for MRP (U, beta, v_beta, mu).
  
    Parameters: 
    df: Merged, filtered, and annotated dataframe containing summary statistics.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    key: Variant/gene name.
    sigma_m_type: One of "sigma_m_var"/"sigma_m_1"/"sigma_m_005". Dictates variant scaling factor by functional annotation.
    R_study: R_study matrix to use for analysis (independent/similar).
    R_phen: R_phen matrix to use for analysis (independent/similar).
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to use for analysis.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    M: Number of variants within the gene block if agg_type is "gene"; if "variant", 1.
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies and phenotypes. Used to calculate v_beta.
  
    Returns: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K) dictating correlation structures, adjusted for missingness.
    beta: A S*M*K x 1 vector of effect sizes.
    v_beta: A (S*M*K x S*M*K) matrix of variances of effect sizes.
    mu: A mean of genetic effects, size of beta (NOTE: default is 0, can change in the code below).
  
    """

    subset_df = (
        df[df["gene_symbol"] == key] if agg_type == "gene" else df[df["V"] == key]
    )
    sigma_m = subset_df[sigma_m_type].tolist()
    diag_sigma_m = np.diag(np.atleast_1d(np.array(sigma_m)))
    R_var = np.diag(np.ones(M)) if R_var_model == "independent" else np.ones((M, M))
    S_var = np.dot(np.dot(diag_sigma_m, R_var), diag_sigma_m)
    beta_list, se_list = generate_beta_se(subset_df, pops, phenos)
    beta = np.array(beta_list).reshape(-1, 1)
    se = np.array(se_list)
    omega = np.kron(err_corr, np.diag(np.ones(M)))
    U = np.kron(np.kron(R_study, R_phen), S_var)
    U, omega, beta, se = adjust_for_missingness(U, omega, beta, se, beta_list)
    diag_se = np.diag(se)
    v_beta = np.dot(np.dot(diag_se, omega), diag_se)
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
    err_corr,
    R_var_model,
    analysis,
    sigma_m_type,
    conserved_columns,
    agg_type,
):

    """ 
    Runs MRP with the given parameters.
  
    Parameters: 
    dfs: List of dataframes that have been filtered and annotated for analysis; need to be merged.
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    R_study: R_study matrix to use for analysis (independent/similar).
    R_study_model: String ("independent"/"similar") corresponding to R_study.
    R_phen: R_phen matrix to use for analysis (independent/similar).
    R_phen_model: String ("independent"/"similar") corresponding to R_phen.
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies and phenotypes. Used to calculate v_beta.
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to use for analysis.
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.
    sigma_m_type: One of "sigma_m_var"/"sigma_m_1"/"sigma_m_005". Dictates variant scaling factor by functional annotation.
    conserved_columns: Columns in every annotated summary statistic file; basis of merges.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.

    Returns: 
    bf_df: Dataframe with two columns: agg_type and log_10 Bayes Factor. 
  
    """

    outer_merge = partial(pd.merge, on=conserved_columns, how="outer")
    df = reduce(outer_merge, dfs)
    m_dict = (
        df.groupby("gene_symbol").size()
        if agg_type == "gene"
        else df.groupby("V").size()
    )
    bf_dict = {}
    for i, (key, value) in enumerate(m_dict.items()):
        if i % 1000 == 0:
            print("Done " + str(i) + " " + agg_type + "s out of " + str(len(m_dict)))
        M = value
        U, beta, v_beta, mu = calculate_all_params(
            df,
            pops,
            phenos,
            key,
            sigma_m_type,
            R_study,
            R_phen,
            R_var_model,
            agg_type,
            M,
            err_corr,
        )
        bf = return_BF(U, beta, v_beta, mu, key, agg_type)
        bf_dict[key] = bf
    bf_df = (
        pd.DataFrame.from_dict(bf_dict, orient="index")
        .reset_index()
        .rename(
            columns={
                "index": agg_type,
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

    """ 
    Filters a set of dataframes that have been read in based on functional consequence.
  
    Dependent on the variant filter that is dictated by the analysis.
  
    Parameters: 
    sumstats_files: The list of summary statistic files that have been read in and annotated.
    variant_filter: The variant filter dictated by the analysis ("ptv"/"pav"/"pcv").
  
    Returns: 
    analysis_files: A list of summary statistic files ready for MRP.
  
    """

    analysis_files = []
    for df in sumstats_files:
        if variant_filter == "ptv":
            df = df[df.category == "ptv"]
        elif variant_filter == "pav":
            df = df[(df.category == "ptv") | (df.category == "pav")]
        analysis_files.append(df)
    return analysis_files


def rename_columns(df, conserved_columns, pop, pheno):

    """ 
    Renames columns such that information on population/study and phenotype is available in the resultant df. 
  
    Additionally checks if the header contains "LOG(OR)_SE" instead of "SE".
  
    Parameters: 
    df: Input dataframe (from summary statistics).
    conserved_columns: These columns are not renamed so that we can merge summary statistics across phenotypes and studies.
    pop: The study from which the current summary statistic df comes from
    pheno: The phenotype from which the current summary statistic df comes from.
  
    Returns: 
    df: A df with adjusted column names, e.g., "OR_white_british_cancer1085".
  
    """

    if "LOG(OR)_SE" in df.columns:
        df.rename(columns={"LOG(OR)_SE": "SE"}, inplace=True)
    columns_to_rename = list(set(df.columns) - set(conserved_columns))
    renamed_columns = [(x + "_" + pop + "_" + pheno) for x in columns_to_rename]
    df.rename(columns=dict(zip(columns_to_rename, renamed_columns)), inplace=True)
    return df


def get_beta(df, pop, pheno):

    """
    Retrieves betas from one (pop, pheno) tuple using non-significant, non-missing variants.
  
    Parameters:
    df: Merged dataframe containing summary statistics.
    pop: Population of interest.
    pheno: Phenotype of interest.
  
    Returns: 
    beta: List of effect sizes from the specified (pop, pheno) tuple; used to compute correlation.
  
    """

    if "BETA" + "_" + pop + "_" + pheno in df.columns:
        return list(df["BETA" + "_" + pop + "_" + pheno])
    elif "OR" + "_" + pop + "_" + pheno in df.columns:
        return np.log(list(df["OR" + "_" + pop + "_" + pheno]))
    else:
        return None


def get_betas(df, pop1, pheno1, pop2, pheno2):

    """
    Retrieves betas from a pair of (pop, pheno) tuples using non-significant, non-missing variants.
  
    Parameters: 
    df: Merged dataframe containing summary statistics.
    pop1: First population.
    pheno1: First phenotype.
    pop2: Second population.
    pheno2: Second phenotype.
  
    Returns: 
    beta1: List of effect sizes from the first (pop, pheno) tuple; used to compute correlation.
    beta2: List of effect sizes from the second (pop, pheno) tuple; used to compute correlation.
  
    """

    df = df[
        (df["P_" + pop1 + "_" + pheno1] >= 0.01)
        & (df["P_" + pop2 + "_" + pheno2] >= 0.01)
    ]
    beta1 = get_beta(df, pop1, pheno1)
    beta2 = get_beta(df, pop2, pheno2)
    return beta1, beta2


def build_err_corr(dfs, pops, phenos, S, K):

    """ 
    Builds a matrix of correlations of errors across studies and phenotypes.
  
    Parameters: 
    dfs: List of dataframes that contain summary statistics.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    S: Number of populations/studies.
    K: Number of GBE phenotypes.

    Returns:
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies and phenotypes. Used to calculate v_beta.
  
    """

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
    ]

    outer_merge = partial(pd.merge, on=conserved_columns, how="outer")
    df = reduce(outer_merge, dfs)
    null_variants = [
        "regulatory_region_variant",
        "intron_variant",
        "intergenic_variant",
        "downstream_gene_variant",
        "mature_miRNA_variant",
        "non_coding_transcript_exon_variant",
        "upstream_gene_variant",
        "NA",
        "NMD_transcript_variant",
        "synonymous_variant",
    ]
    # Get only null variants
    df = df[df.most_severe_consequence.isin(null_variants)]
    # Get only common variants
    df = df[df.wb_maf >= 0.01]

    # Filter out any rows that have NA values
    df = df[
        conserved_columns
        + [
            col
            for col in df.columns
            if (("BETA" in col) or ("OR" in col) or ("P" in col))
        ]
    ]
    df = df.dropna()

    err_corr = np.zeros((S * K, S * K))
    for i, pop1 in enumerate(pops):
        for j, pheno1 in enumerate(phenos):
            for x, pop2 in enumerate(pops):
                for y, pheno2 in enumerate(phenos):
                    # Location in matrix
                    a = K * i + j
                    b = K * x + y
                    # If in lower triangle, do not compute; symmetric matrix
                    if a > b:
                        err_corr[a, b] = err_corr[b, a]
                    else:
                        beta1, beta2 = get_betas(df, pop1, pheno1, pop2, pheno2)
                        err_corr[a, b] = pearsonr(beta1, beta2)[0]
    return err_corr


def collect_and_filter(pops, phenos, datasets, conserved_columns):

    """ 
    Collects the summary statistics of interest and applies filters. 
  
    Reads in summary statistics; 
    Joins on metadata files that contain gene symbol, MAF, and consequence;
    Sets sigma values associated with different types of consequences (PTV/PAV/PCV).
    
    Parameters: 
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    datasets: Unique list of datasets ("cal"/"exome") to use for analysis.
    conserved_columns: Columns in every annotated summary statistic file; basis of merges.
  
    Returns: 
    filtered_sumstats_files: Annotated and filtered summary statistic files.
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies and phenotypes. Used to calculate v_beta.
  
    """

    print("Reading in summary stats for:")
    print("")
    print("Populations: " + ", ".join(pops))
    print("Phenotypes: " + ", ".join(phenos))
    print("Datasets: " + ", ".join(datasets))
    print("")

    file_paths, sumstat_files = read_in_summary_stats(
        pops, phenos, datasets, conserved_columns
    )

    metadata_sumstat_files = []
    exome_metadata = pd.read_csv(
        "/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-gene_consequence_wb_maf_final.tsv",
        sep="\t",
    )
    cal_metadata = pd.read_csv(
        "/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-gene_consequence_wb_maf_final.tsv",
        sep="\t",
    )

    for file_path, sumstat_file in zip(file_paths, sumstat_files):
        # Merge metadata
        if "exome" in file_path:
            sumstat_file = sumstat_file.merge(exome_metadata)
        elif "cal" in file_path:
            sumstat_file = sumstat_file.merge(cal_metadata)
        metadata_sumstat_files.append(sumstat_file)

    # Sample common variants, stuff in filter + synonymous
    print("")
    print("Building matrix of correlations of errors...")
    err_corr = build_err_corr(
        metadata_sumstat_files, pops, phenos, len(pops), len(phenos)
    )
    filtered_sumstat_files = []

    print("")
    print("Filtering sumstats on MAF, setting sigmas, and filtering on consequence...")
    print("")

    for sumstat_file in metadata_sumstat_files:
        filtered_sumstat_file = sumstat_file[
            (sumstat_file.wb_maf <= 0.01) & (sumstat_file.wb_maf > 0)
        ]
        filtered_sumstat_file = set_sigmas(filtered_sumstat_file)
        filtered_sumstat_files.append(filtered_sumstat_file)
    return filtered_sumstat_files, err_corr


def return_input_args(args):

    """ 
    Further parses the command-line input.
  
    Makes all lists unique; calculates S and K; and creates lists of appropriate matrices.
  
    Parameters: 
    args: Command-line arguments that have been parsed by the parser.
  
    Returns: 
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    R_study_list: Unique list of R_study matrices to use for analysis.
    R_study_models: Unique strings ("independent"/"similar") corresponding to each matrix in R_study_list.
    R_phen_list: Unique list of R_phen matrices to use for analysis.
    R_phen_models: Unique strings ("independent"/"similar") corresponding to each matrix in R_phen_list.
    R_var_models: Unique strings ("independent"/"similar") corresponding to R_var matrices to use for analysis.
    agg: Unique list of aggregation units ("gene"/"variant") to use for analysis.
    sigma_m_types: Unique list of sigma_m types ("sigma_m_var"/"sigma_m_1"/"sigma_m_005") to use for analysis.
    variant_filters: Unique list of variant filters ("ptv"/"pav"/"pcv") to use for analysis.
    datasets: Unique list of datasets ("cal"/"exome") to use for analysis.
  
    """

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


def print_banner():

    """ 
    Prints ASCII Art Banner + Author Info.
  
    """

    cprint(" __  __ ____  ____", "red")
    cprint("|  \/  |  _ \|  _ \\", "red")
    cprint("| |\/| | |_) | |_) |", "red")
    cprint("| |  | |  _ <|  __/ ", "red")
    cprint("|_|  |_|_| \_\_|  ", "red")
    print("")
    cprint("Production Author:", "green")
    print("Guhan Ram Venkataraman, B.S.H.")
    print("Ph.D. Candidate | Biomedical Informatics")
    print("")
    print("Contact: guhan@stanford.edu")
    print(
        "URL: https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_production.py"
    )
    print("")
    print("Methods Developers:", "green")
    print("Manuel A. Rivas, Ph.D.; Matti Pirinen, Ph.D.")
    print("Rivas Lab | Stanford University")
    print("")


def print_params(
    analysis, R_study_model, R_phen_model, R_var_model, agg_type, sigma_m_type
):

    """ 
    Provides a text overview of each analysis in the terminal.
  
    Parameters: 
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.
    R_study_model: One of "independent"/"similar". Dictates correlation structure across studies.
    R_phen_model: One of "independent"/"similar". Dictates correlation structure across phenotypes.
    R_var_model: One of "independent"/"similar". Dictates correlation structure across variants.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    sigma_m_type: One of "sigma_m_var"/"sigma_m_1"/"sigma_m_005". Dictates variant scaling factor by functional annotation.
  
    """

    print("Analysis: " + analysis)
    print("R_study model: " + R_study_model)
    print("R_phen model: " + R_phen_model)
    print("R_var model: " + R_var_model)
    print("Aggregation by: " + agg_type)
    print("Variant weighting factor: " + sigma_m_type)
    print("")


def initialize_parser(valid_phenos):

    """
    Parses inputs using argparse. 
    
    Parameters: 
    valid_phenos: List of valid GBE phenotypes read in from phenotype_info.tsv.

    """
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
        default=["similar"],
        dest="R_study_models",
        help="type of model across studies. options: independent, similar (default: similar). can run both.",
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
        choices=["sigma_m_var", "sigma_m_1", "sigma_m_005"],
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
    return parser


def output_file(bf_dfs, agg_type, dataset, pops, phenos):

    """ 
    Outputs a file containing aggregation unit and Bayes Factors. 
    
    Parameters: 
    bf_dfs: List of dataframes containing Bayes Factors from each anOne of Unique set of GBE phenotypes to use for analysis.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    dataset: One of "cal"/"exome".
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.

    """

    outer_merge = partial(pd.merge, on=agg_type, how="outer")
    out_df = reduce(outer_merge, bf_dfs)
    out_file = (
        "/oak/stanford/groups/mrivas/ukbb24983/"
        + dataset
        + "/mrp/"
        + "_".join(pops)
        + "_"
        + "_".join(phenos)
        + "_"
        + agg_type
        + ".tsv"
    )
    out_df.to_csv(out_file, sep="\t", index=False)
    print("")
    print("Results written to " + out_file + ".")


def loop_through_parameters(
    datasets,
    agg,
    variant_filters,
    S,
    R_study_list,
    R_study_models,
    pops,
    K,
    R_phen_list,
    R_phen_models,
    phenos,
    R_var_models,
    sigma_m_types,
    sumstats_files,
    err_corr,
    conserved_columns,
):

    """ 
    Loops through parameters specified through command line (or defaults). 

    Parameters: 
    datasets: Unique list of datasets ("cal"/"exome") to use for analysis.
    agg: Unique list of aggregation units ("gene"/"variant") to use for analysis.
    variant_filters: Unique list of variant filters ("ptv"/"pav"/"pcv") to use for analysis.
    S: Number of populations/studies.
    R_study_list: Unique list of R_study matrices to use for analysis.
    R_study_models: Unique strings ("independent"/"similar") corresponding to each matrix in R_study_list.
    pops: Unique set of populations (studies) to use for analysis.
    K: Number of GBE phenotypes.
    R_phen_list: Unique list of R_phen matrices to use for analysis.
    R_phen_models: Unique strings ("independent"/"similar") corresponding to each matrix in R_phen_list.
    phenos: Unique set of GBE phenotypes to use for analysis.
    R_var_models: Unique strings ("independent"/"similar") corresponding to R_var matrices to use for analysis.
    sigma_m_types: Unique list of sigma_m types ("sigma_m_var"/"sigma_m_1"/"sigma_m_005") to use for analysis.
    sumstats_files: List of dataframes containing filtered summary statistic files.
    err_corr: Matrix of correlation of errors across studies and phenotypes.
    conserved_columns: Columns in every annotated summary statistic file; basis of merges.
  
    """

    for dataset in datasets:
        for agg_type in agg:
            bf_dfs = []
            for analysis in variant_filters:
                analysis_files = filter_category(sumstat_files, analysis)
                for R_study, R_study_model in zip(R_study_list, R_study_models):
                    for R_phen, R_phen_model in zip(R_phen_list, R_phen_models):
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
                                    err_corr,
                                    R_var_model,
                                    analysis,
                                    sigma_m_type,
                                    conserved_columns,
                                    agg_type,
                                )
                                bf_df = bf_df.sort_values(
                                    by=list(set(bf_df.columns) - set([agg_type])),
                                    ascending=False,
                                )
                                bf_dfs.append(bf_df)
            output_file(bf_dfs, agg_type, dataset, pops, phenos)


if __name__ == "__main__":

    """ 
    Runs MRP analysis on GBE summary statistics with the parameters specified by the command line.

    """

    with open("../05_gbe/phenotype_info.tsv", "r") as phe_file:
        valid_phenos = [line.split()[0] for line in phe_file][1:]

    parser = initialize_parser(valid_phenos)
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
    from scipy.stats.stats import pearsonr
    from collections import defaultdict
    import subprocess
    import math
    import os
    from termcolor import cprint

    print_banner()

    S, K, pops, phenos, R_study_list, R_study_models, R_phen_list, R_phen_models, R_var_models, agg, sigma_m_types, variant_filters, datasets = return_input_args(
        args
    )

    # These columns are in every annotated summary statistic file and will be the basis of merges
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

    sumstat_files, err_corr = collect_and_filter(
        pops, phenos, datasets, conserved_columns
    )

    loop_through_parameters(
        datasets,
        agg,
        variant_filters,
        S,
        R_study_list,
        R_study_models,
        pops,
        K,
        R_phen_list,
        R_phen_models,
        phenos,
        R_var_models,
        sigma_m_types,
        sumstats_files,
        err_corr,
        conserved_columns,
    )
