from __future__ import division
import argparse


def is_pos_def_and_full_rank(X):

    """ 
    Ensures a matrix is positive definite and full rank.
  
    Keep diagonals and multiples every other cell by .99.
  
    Parameters: 
    X: Matrix to verify.
  
    Returns: 
    X: Verified (and, if applicable, adjusted) matrix.
  
    """

    i = 0
    X = np.matrix(X)
    if (np.all(np.linalg.eigvals(X) > 0)) and (np.linalg.matrix_rank(X) == len(X)):
        return X
    else:
        while not (
            (np.all(np.linalg.eigvals(X) > 0)) and (np.linalg.matrix_rank(X) == len(X))
        ):
            X = 0.99 * X + 0.01 * np.diag(np.diag(X))
        return X


def safe_inv(X, matrix_name, block, agg_type):

    """ 
    Safely inverts a matrix, or returns NaN.
  
    Parameters: 
    X: Matrix to invert.
    matrix_name: One of "U"/"v_beta" - used to print messages when inversion fails.
    block: Name of the aggregation block (gene/variant). 
        Used to print messages when inversion fails.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
        Used to print messages when inversion fails.
  
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


def farebrother(quad_T, d, fb):

    """ 
    Farebrother method from CompQuadForm.
  
    Parameters: 

    quad_T: Value point at which distribution function is to be evaluated.
    d: Distinct non-zero characteristic root(s) of A*Sigma. 
    fb: Farebrother R method (rpy2 object).
  
    Returns: 
    p_value: Farebrother p-value.
  
    """

    res = fb(quad_T, d)
    return np.asarray(res)[0]


def davies(quad_T, d, dm):

    """ 
    Davies method from CompQuadForm.
  
    Parameters: 

    quad_T: Value point at which distribution function is to be evaluated.
    d: Distinct non-zero characteristic root(s) of A*Sigma. 
    dm: Davies R method (rpy2 object).
  
    Returns: 
    p_value: Davies p-value.
  
    """

    res = dm(quad_T, d)
    return np.asarray(res)[0]


def imhof(quad_T, d, im):

    """ 
    Imhof method from CompQuadForm.
  
    Parameters: 

    quad_T: Value point at which distribution function is to be evaluated.
    d: Distinct non-zero characteristic root(s) of A*Sigma. 
    im: Imhof R method (rpy2 object).
  
    Returns: 
    p_value: Imhof p-value.
  
    """

    res = im(quad_T, d)
    return np.asarray(res)[0]


def initialize_r_objects():

    """ 
    Initializes Farebrother, Davies, and Imhof R methods as rpy2 objects.
  
    Returns: 
    fb: Farebrother R method (rpy2 object).
    dm: Davies R method (rpy2 object).
    im: Imhof R method (rpy2 object).
  
    """

    robjects.r(
        """
    require(MASS)
    require(CompQuadForm)
    farebrother.method <- function(quadT, d, h = rep(1, length(d)), delta = rep(0, length(d)), maxiter = 100000, epsilon = 10^-16, type = 1) {
        return(farebrother(quadT, d, h, delta, maxit = as.numeric(maxiter), eps = as.numeric(epsilon), mode = as.numeric(type))$Qq)
    }
    """
    )
    robjects.r(
        """
    require(MASS)
    require(CompQuadForm)
    imhof.method <- function(quadT, d, h = rep(1, length(d)), delta = rep(0, length(d)), epsilon = 10^-16, lim = 10000) {
        return(imhof(quadT, d, h, delta, epsabs = as.numeric(epsilon), epsrel = as.numeric(epsilon), limit=as.numeric(lim))$Qq)
    }
    """
    )
    robjects.r(
        """
    require(MASS)
    require(CompQuadForm)
    davies.method <- function(quadT, d, h = rep(1, length(d)), delta = rep(0, length(d)), sig = 0, limit = 10000, accuracy = 0.0001) {
        return(davies(quadT, d, h, delta, sigma=as.numeric(sig), lim = as.numeric(limit), acc = as.numeric(accuracy))$Qq)
    }
    """
    )
    fb = robjects.globalenv["farebrother.method"]
    fb = robjects.r["farebrother.method"]
    im = robjects.globalenv["imhof.method"]
    im = robjects.r["imhof.method"]
    dm = robjects.globalenv["davies.method"]
    dm = robjects.r["davies.method"]
    return fb, dm, im


def return_BF_pvals(beta, U, v_beta, v_beta_inv, fb, dm, im, methods):

    """ 
    Computes a p-value from the quadratic form that is subsumed by the Bayes Factor.
  
    Parameters: 

    beta: Effect size vector without missing data.
    U: Kronecker product of the three matrices (S*M*K x S*M*K)
        dictating correlation structures; no missing data.
    v_beta: Diagonal matrix of variances of effect sizes without missing data.
    v_beta_inv: Inverse of v_beta.
    fb: Farebrother R method (rpy2 object).
    dm: Davies R method (rpy2 object).
    im: Imhof R method (rpy2 object).
    methods: List of method(s) to apply to our data.
  
    Returns: 
    p_values: List of p-values corresponding to each method specified as input.
  
    """

    n = beta.shape[0]
    A = v_beta + U
    A = is_pos_def_and_full_rank(A)
    A_inv = np.linalg.inv(A)
    quad_T = beta.T * (v_beta_inv - A_inv) * beta
    B = is_pos_def_and_full_rank(np.eye(n) - A_inv * v_beta)
    d = np.linalg.eig(B)[0]
    d = [i for i in d if i > 0.01]
    methods = ["imhof", "davies", "farebrother"]
    p_values = []
    for method in methods:
        if method == "farebrother":
            p_value = farebrother(quad_T, d, fb)
        elif method == "davies":
            p_value = davies(quad_T, d, dm)
        elif method == "imhof":
            p_value = imhof(quad_T, d, im)
        p_value = max(0, min(1, p_value))
        p_values.append(p_value)
    return p_values


def compute_posterior_probs(log10BF, prior_odds_list):

    """ 
    Computes posterior probability given prior odds and a log10 Bayes Factor.
  
    Parameters: 
    log10BF: log10 Bayes Factor of given association.
    prior_odds_list: List of assumed prior odds.
  
    Returns: 
    posterior_prob: The posterior probability of the event
        given the prior odds and the Bayes Factor.
  
    """

    BF = 10 ** (log10BF)
    posterior_odds_list = [prior_odds * BF for prior_odds in prior_odds_list]
    posterior_probs = [
        (posterior_odds / (1 + posterior_odds))
        for posterior_odds in posterior_odds_list
    ]
    return posterior_probs


def return_BF(
    U, beta, v_beta, mu, block, agg_type, prior_odds_list, p_value_methods, fb, dm, im
):

    """ 
    Given quantities calculated previously and the inputs, returns the associated 
        Bayes Factor.
  
    Parameters: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K)
        dictating correlation structures; no missing data.
    beta: Effect size vector without missing data.
    v_beta: Diagonal matrix of variances of effect sizes without missing data.
    mu: A mean of genetic effects, size of beta 
        (NOTE: default is 0, can change in the code below).
    block: Name of the aggregation block (gene/variant).
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    prior_odds_list: List of prior odds used as assumptions to calculate 
        posterior probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from 
        Bayes Factors.
    fb, dm, im: initialized R functions for Farebrother, Davies, and Imhof methods. 
        NoneType if p_value_methods is [].
  
    Returns:, [] 
    log10BF: log_10 Bayes Factor (ratio of marginal likelihoods of alternative model,
        which accounts for priors, and null).
    posterior_probs: List of posterior probabilities corresponding to each prior odds
         in prior_odds_list.
    p_values: List of p-values corresponding to each method in p_value_methods.
  
    """

    v_beta = is_pos_def_and_full_rank(v_beta)
    v_beta_inv = safe_inv(v_beta, "v_beta", block, agg_type)
    U = is_pos_def_and_full_rank(U)
    U_inv = safe_inv(U, "U", block, agg_type)
    if v_beta_inv is not np.nan and U_inv is not np.nan:
        A2 = U_inv + v_beta_inv
        b2 = v_beta_inv * beta
        try:
            Abinv = np.linalg.lstsq(A2, b2, rcond=-1)[0]
        except LinAlgError as err:
            return np.nan, [], []
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
        posterior_probs = (
            compute_posterior_probs(log10BF, prior_odds_list) if prior_odds_list else []
        )
        p_values = (
            return_BF_pvals(beta, U, v_beta, v_beta_inv, fb, dm, im, p_value_methods)
            if p_value_methods
            else []
        )
        return log10BF, posterior_probs, p_values
    else:
        return np.nan, [], []


def delete_rows_and_columns(X, indices_to_remove):

    """ 
    Helper function to delete rows and columns from a matrix.
  
    Parameters: 
    X: Matrix that needs adjustment.
    indices_to_remove: Rows and columns to be deleted.
  
    Returns: 
    X: Smaller matrix that has no missing data.
  
    """

    X = np.delete(X, indices_to_remove, axis=0)
    X = np.delete(X, indices_to_remove, axis=1)
    return X


def adjust_for_missingness(U, omega, beta, se, beta_list):

    """ 
    Deletes rows and columns where we do not have effect sizes/standard errors.

    Calls method delete_rows_and_columns, a helper function that calls the numpy
        command.
  
    Parameters: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K) 
        dictating correlation structures; may relate to missing data.
    omega: (S*M*K x S*M*K) matrix that contains correlation of errors
        across variants, studies, and phenotypes. 
    beta: Vector of effect sizes within the unit of aggregation;
        may contain missing data.
    se: Vector of standard errors within the unit of aggregation;
        may contain missing data.
    beta_list: List of effect sizes within the unit of aggregation;
        may contain missing data.
  
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
    subset_df: Slice of the original dataframe that encompasses the current unit of 
        aggregation (gene/variant).
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
            if "BETA" + "_" + pop + "_" + pheno in subset_df.columns:
                beta_list.extend(list(subset_df["BETA" + "_" + pop + "_" + pheno]))
            elif "OR" + "_" + pop + "_" + pheno in subset_df.columns:
                beta_list.extend(list(subset_df["OR" + "_" + pop + "_" + pheno]))
            se_list.extend(list(subset_df["SE" + "_" + pop + "_" + pheno]))
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
    sigma_m_type: One of "sigma_m_var"/"sigma_m_1"/"sigma_m_005".
        Dictates variant scaling factor by functional annotation.
    R_study: R_study matrix to use for analysis (independent/similar).
    R_phen: R_phen matrix to use for analysis (empirically calculated).
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to 
        use for analysis.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    M: Number of variants within the gene block if agg_type is "gene"; if "variant", 1.
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies 
        and phenotypes. Used to calculate v_beta.
  
    Returns: 
    U: Kronecker product of the three matrices (S*M*K x S*M*K)
        dictating correlation structures, adjusted for missingness.
    beta: A S*M*K x 1 vector of effect sizes.
    v_beta: A (S*M*K x S*M*K) matrix of variances of effect sizes.
    mu: A mean of genetic effects, size of beta
        (NOTE: default is 0, can change in the code below).
  
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


def output_file(bf_dfs, agg_type, dataset, pops, phenos, maf_thresh):

    """ 
    Outputs a file containing aggregation unit and Bayes Factors. 
    
    Parameters: 
    bf_dfs: List of dataframes containing Bayes Factors from each analysis.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    dataset: One of "cal"/"exome".
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    maf_thresh: Maximum MAF of variants in this run.

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
        + "_"
        + str(maf_thresh)
        + ".tsv"
    )
    out_df.to_csv(out_file, sep="\t", index=False)
    print("")
    print(Fore.RED + "Results written to " + out_file + "." + Style.RESET_ALL)
    print("")


def run_mrp(
    dfs,
    S,
    K,
    pops,
    phenos,
    R_study,
    R_study_model,
    R_phen,
    err_corr,
    R_var_model,
    analysis,
    sigma_m_type,
    conserved_columns,
    agg_type,
    prior_odds_list,
    p_value_methods,
):

    """ 
    Runs MRP with the given parameters.
  
    Parameters: 
    dfs: List of dataframes that have been filtered and annotated for analysis;
        need to be merged.
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    R_study: R_study matrix to use for analysis (independent/similar).
    R_study_model: String ("independent"/"similar") corresponding to R_study.
    R_phen: R_phen matrix to use for analysis (empirically calculated).
    err_corr: A (S*K x S*K) matrix of correlation of errors across studies and 
        phenotypes. Used to calculate v_beta.
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to 
        use for analysis.
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.
    sigma_m_type: One of "sigma_m_var"/"sigma_m_1"/"sigma_m_005". Dictates variant 
        scaling factor by functional annotation.
    conserved_columns: Columns in every annotated summary statistic file; 
        basis of merges.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    prior_odds_list: List of prior odds used as assumptions to calculate posterior 
        probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from Bayes 
        Factors.

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
    bf_df_columns = [
        agg_type,
        "log_10_BF"
        + "_study_"
        + R_study_model
        + "_var_"
        + R_var_model
        + "_"
        + sigma_m_type
        + "_"
        + analysis,
    ]
    if prior_odds_list:
        bf_df_columns.extend(
            [
                "posterior_prob_w_prior_odds_"
                + str(prior_odds)
                + "_study_"
                + R_study_model
                + "_var_"
                + R_var_model
                + "_"
                + sigma_m_type
                + "_"
                + analysis
                for prior_odds in prior_odds_list
            ]
        )
    if p_value_methods:
        fb, dm, im = initialize_r_objects()
        bf_df_columns.extend(
            [
                "p_value_"
                + p_value_method
                + "_study_"
                + R_study_model
                + "_var_"
                + R_var_model
                + "_"
                + sigma_m_type
                + "_"
                + analysis
                for p_value_method in p_value_methods
            ]
        )
    else:
        fb = dm = im = None
    data = []
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
        bf, posterior_probs, p_values = return_BF(
            U,
            beta,
            v_beta,
            mu,
            key,
            agg_type,
            prior_odds_list,
            p_value_methods,
            fb,
            dm,
            im,
        )
        data.append([key, bf] + posterior_probs + p_values)
    bf_df = pd.DataFrame(data, columns=bf_df_columns)
    return bf_df


def print_params(
    analysis,
    R_study_model,
    R_var_model,
    agg_type,
    sigma_m_type,
    maf_thresh,
    prior_odds_list,
    p_value_methods,
):

    """ 
    Provides a text overview of each analysis in the terminal.
  
    Parameters: 
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.
    R_study_model: One of "independent"/"similar". Dictates correlation structure
        across studies.
    R_var_model: One of "independent"/"similar". Dictates correlation structure
        across variants.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    sigma_m_type: One of "sigma_m_var"/"sigma_m_1"/"sigma_m_005". Dictates variant 
        scaling factor by functional annotation.
    maf_thresh: Maximum MAF of variants in this run.
    prior_odds_list: List of prior odds used as assumptions to calculate posterior 
        probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from Bayes
        Factors.
  
    """

    print("")
    print(Fore.YELLOW + "Analysis: " + Style.RESET_ALL + analysis)
    print(Fore.YELLOW + "R_study model: " + Style.RESET_ALL + R_study_model)
    print(Fore.YELLOW + "R_var model: " + Style.RESET_ALL + R_var_model)
    print(Fore.YELLOW + "Aggregation by: " + Style.RESET_ALL + agg_type)
    print(Fore.YELLOW + "Variant weighting factor: " + Style.RESET_ALL + sigma_m_type)
    print(Fore.YELLOW + "MAF threshold: " + Style.RESET_ALL + str(maf_thresh))
    if prior_odds_list:
        print(
            Fore.YELLOW
            + "Prior odds: "
            + Style.RESET_ALL
            + ", ".join([str(prior_odd) for prior_odd in prior_odds_list])
        )
    if p_value_methods:
        print(
            Fore.YELLOW
            + "Methods for p-value generation: "
            + Style.RESET_ALL
            + ", ".join(p_value_methods)
        )
    print("")


def filter_category(sumstats_files, variant_filter):

    """ 
    Filters a set of dataframes that have been read in based on functional consequence.
  
    Dependent on the variant filter that is dictated by the analysis.
  
    Parameters: 
    sumstats_files: The list of summary statistic files that have been read
        in and annotated.
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


def loop_through_parameters(
    dataset,
    agg,
    variant_filters,
    S,
    R_study_list,
    R_study_models,
    pops,
    K,
    R_phen,
    phenos,
    R_var_models,
    sigma_m_types,
    sumstat_files,
    err_corr,
    conserved_columns,
    maf_thresh,
    prior_odds_list,
    p_value_methods,
):

    """ 
    Loops through parameters specified through command line (or defaults). 

    Parameters: 
    dataset: One of ("cal"/"exome") to use for analysis.
    agg: Unique list of aggregation units ("gene"/"variant") to use for analysis.
    variant_filters: Unique list of variant filters ("ptv"/"pav"/"pcv") to use 
        for analysis.
    S: Number of populations/studies.
    R_study_list: Unique list of R_study matrices to use for analysis.
    R_study_models: Unique strings ("independent"/"similar") corresponding to each 
        matrix in R_study_list.
    pops: Unique set of populations (studies) to use for analysis.
    K: Number of GBE phenotypes.
    R_phen: R_phen matrix to use for analysis (empirically calculated).
    phenos: Unique set of GBE phenotypes to use for analysis.
    R_var_models: Unique strings ("independent"/"similar") corresponding to R_var 
        matrices to use for analysis.
    sigma_m_types: Unique list of sigma_m types ("sigma_m_var"/"sigma_m_1"/"sigma_m_005")
        to use for analysis.
    sumstats_files: List of dataframes containing filtered summary statistic files.
    err_corr: Matrix of correlation of errors across studies and phenotypes.
    conserved_columns: Columns in every annotated summary statistic file; 
        basis of merges.
    maf_thresh: Maximum MAF of variants in this run.
    prior_odds_list: List of prior odds used as assumptions to calculate posterior 
        probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from Bayes 
        Factors.
  
    """

    print(
        Fore.YELLOW
        + "Running MRP across parameters for maf_thresh "
        + str(maf_thresh)
        + "..."
        + Style.RESET_ALL
    )
    for agg_type in agg:
        bf_dfs = []
        # If not aggregating, then R_var choice does not affect BF (just a 1x1 matrix, [1])
        if (agg_type == "variant") and (len(R_var_models) > 1):
            print(Fore.YELLOW)
            print("Since we are not aggregating, R_var is just [1].")
            print(Style.RESET_ALL)
            R_var_models = ["independent"]
        for analysis in variant_filters:
            analysis_files = filter_category(sumstat_files, analysis)
            for R_study, R_study_model in zip(R_study_list, R_study_models):
                for R_var_model in R_var_models:
                    for sigma_m_type in sigma_m_types:
                        print_params(
                            analysis,
                            R_study_model,
                            R_var_model,
                            agg_type,
                            sigma_m_type,
                            maf_thresh,
                            prior_odds_list,
                            p_value_methods,
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
                            err_corr,
                            R_var_model,
                            analysis,
                            sigma_m_type,
                            conserved_columns,
                            agg_type,
                            prior_odds_list,
                            p_value_methods,
                        )
                        bf_df = bf_df.sort_values(
                            by=list(set(bf_df.columns) - set([agg_type])),
                            ascending=False,
                        )
                        bf_dfs.append(bf_df)
        output_file(bf_dfs, agg_type, dataset, pops, phenos, maf_thresh)


def set_sigmas(df):

    """ 
    Assigns appropriate sigmas to appropriate variants by annotation.
  
    Filters out variants not of interest;
    Sets sigmas by functional annotation;
    Additionally adds two extra columns for two other standard choices of a uniform 
        sigma (1 and 0.05).
  
    Parameters: 
    df: Merged dataframe containing all variants across all studies and phenotypes.
  
    Returns: 
    df: Merged dataframe with three additional columns:
        sigma_m_var: Column of sigma values (mapped to functional annotation via the 
            lists inside this method).
            NOTE: One can change the sigmas associated with each type of variant by 
                adjusting the values within this method.
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
    sigma_m_pc = 0.03
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


def filter_and_set_sigmas(metadata_sumstat_files, maf_thresh):

    """
    Filters each dataframe in the list provided given a MAF threshold.

    Parameters:
    metadata_sumstat_files: List of dataframes that have been joined with metadata 
        containing MAF.
    maf_thresh: Upper threshold at which MAF will be cut off.

    Returns:
    filtered_sumstat_files: List of dataframes with only variants with MAF below the 
        threshold specified.

    """

    filtered_sumstat_files = []

    print("")
    print(
        Fore.BLUE
        + "Filtering sumstats on MAF, setting sigmas, and filtering on consequence..."
        + Style.RESET_ALL
    )
    print("")

    for sumstat_file in metadata_sumstat_files:
        filtered_sumstat_file = sumstat_file[
            (sumstat_file.wb_maf <= maf_thresh) & (sumstat_file.wb_maf > 0)
        ]
        filtered_sumstat_file = set_sigmas(filtered_sumstat_file)
        filtered_sumstat_files.append(filtered_sumstat_file)
    return filtered_sumstat_files


def get_beta(df, pop, pheno):

    """
    Retrieves betas from one (pop, pheno) tuple using non-significant, 
        non-missing variants.
  
    Parameters:
    df: Merged dataframe containing summary statistics.
    pop: Population of interest.
    pheno: Phenotype of interest.
  
    Returns: 
    beta: List of effect sizes from the specified (pop, pheno) tuple; 
        used to compute correlation.
  
    """

    if "BETA" + "_" + pop + "_" + pheno in df.columns:
        return list(df["BETA" + "_" + pop + "_" + pheno])
    elif "OR" + "_" + pop + "_" + pheno in df.columns:
        return np.log(list(df["OR" + "_" + pop + "_" + pheno]))
    return None


def get_betas(df, pop1, pheno1, pop2, pheno2, mode):

    """
    Retrieves betas from a pair of (pop, pheno) tuples using non-significant, 
        non-missing variants.
  
    Parameters: 
    df: Merged dataframe containing summary statistics.
    pop1: First population.
    pheno1: First phenotype.
    pop2: Second population.
    pheno2: Second phenotype.
    mode: One of "null", "sig". Determines whether we want to sample from null or 
        significant variants. Useful for building out correlations of errors and 
        phenotypes respectively.
  
    Returns: 
    beta1: List of effect sizes from the first (pop, pheno) tuple; used to compute 
        correlation.
    beta2: List of effect sizes from the second (pop, pheno) tuple; used to compute
        correlation.
  
    """

    if ("P_" + pop1 + "_" + pheno1 not in df.columns) or (
        "P_" + pop2 + "_" + pheno2 not in df.columns
    ):
        return None, None
    if mode == "null":
        df = df[
            (df["P_" + pop1 + "_" + pheno1] >= 1e-2)
            & (df["P_" + pop2 + "_" + pheno2] >= 1e-2)
        ]
    elif mode == "sig":
        df = df[
            (df["P_" + pop1 + "_" + pheno1] < 1e-7)
            | (df["P_" + pop2 + "_" + pheno2] < 1e-7)
        ]
    beta1 = get_beta(df, pop1, pheno1)
    beta2 = get_beta(df, pop2, pheno2)
    return beta1, beta2


def build_phen_corr(S, K, pops, phenos, df):

    """
    Builds out a matrix of correlations between all phenotypes and studies using:
        - significant (P < 1e-7)
        - common (MAF >= 0.01)
        - LD independent
    SNPs.

    Parameters:
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.

    Returns:
    phen_corr: (S*K x S*K) matrix of correlations between all phenotypes and studies 
        for significant variants. Used to calculate R_phen.

    """

    phen_corr = np.zeros((S * K, S * K))
    for i, pop1 in enumerate(pops):
        for j, pheno1 in enumerate(phenos):
            for x, pop2 in enumerate(pops):
                for y, pheno2 in enumerate(phenos):
                    # Location in matrix
                    a, b = K * i + j, K * x + y
                    # If in lower triangle, do not compute; symmetric matrix
                    if a > b:
                        phen_corr[a, b] = np.nan
                    elif a == b:
                        phen_corr[a, b] = np.nan
                    else:
                        phen_beta1, phen_beta2 = get_betas(
                            df, pop1, pheno1, pop2, pheno2, "sig"
                        )
                        phen_corr[a, b] = (
                            pearsonr(phen_beta1, phen_beta2)[0]
                            if phen_beta1 is not None
                            else np.nan
                        )
    return phen_corr


def build_R_phen(S, K, pops, phenos, df):

    """
    Builds R_phen using phen_corr (calculated using the method directly above this).

    Parameters:
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.

    Returns:
    R_phen: Empirical estimates of genetic correlation across phenotypes.

    """

    if K == 1:
        return np.ones((K, K))
    phen_corr = build_phen_corr(S, K, pops, phenos, df)
    R_phen = np.zeros((K, K))
    for k1 in range(K):
        for k2 in range(K):
            if k1 == k2:
                R_phen[k1, k2] = 1
            elif k1 > k2:
                R_phen[k1, k2] = R_phen[k2, k1]
            else:
                phenos_to_remove = list(set(range(K)) - set([k1, k2]))
                indices_to_remove = []
                for pheno_to_remove in phenos_to_remove:
                    indices_to_remove.extend(
                        range(S * pheno_to_remove, S * (pheno_to_remove + 1))
                    )
                pairwise_corrs = delete_rows_and_columns(phen_corr, indices_to_remove)
                R_phen[k1, k2] = np.nanmedian(pairwise_corrs)
    return R_phen


def build_err_corr(S, K, pops, phenos, df):

    """
    Builds out a matrix of correlations between all phenotypes and studies using:
        - null (i.e. synonymous or functionally uninteresting)
        - not significant (P >= 1e-2)
        - common (MAF >= 0.01)
        - LD independent
    SNPs.

    Parameters:
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.

    Returns:
    err_corr: (S*K x S*K) matrix of correlation of errors across studies and phenotypes 
        for null variants. Used to calculate v_beta.

    """

    err_corr = np.zeros((S * K, S * K))
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
    # Get only null variants to build err_corr
    err_df = df[df.most_severe_consequence.isin(null_variants)]
    for i, pop1 in enumerate(pops):
        for j, pheno1 in enumerate(phenos):
            for x, pop2 in enumerate(pops):
                for y, pheno2 in enumerate(phenos):
                    # Location in matrix
                    a, b = K * i + j, K * x + y
                    # If in lower triangle, do not compute; symmetric matrix
                    if a > b:
                        err_corr[a, b] = err_corr[b, a]
                    elif a == b:
                        err_corr[a, b] = 1
                    else:
                        err_beta1, err_beta2 = get_betas(
                            err_df, pop1, pheno1, pop2, pheno2, "null"
                        )
                        err_corr[a, b] = (
                            pearsonr(err_beta1, err_beta2)[0]
                            if err_beta1 is not None
                            else 0
                        )
    return err_corr


def return_err_and_R_phen(dfs, pops, phenos, S, K):

    """ 
    Builds a matrix of correlations of errors across studies and phenotypes,
        and correlations of phenotypes.
  
    Parameters: 
    dfs: List of dataframes that contain summary statistics.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    S: Number of populations/studies.
    K: Number of GBE phenotypes.

    Returns:
    err_corr: (S*K x S*K) matrix of correlation of errors across studies and phenotypes
        for null variants. Used to calculate v_beta.
    R_phen: Empirical estimates of genetic correlation across phenotypes.
  
    """

    # Sample common variants, stuff in filter + synonymous
    print("")
    print(
        Fore.MAGENTA + "Building matrix of correlations of errors..." + Style.RESET_ALL
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
        "ld_indep",
    ]
    outer_merge = partial(pd.merge, on=conserved_columns, how="outer")
    df = reduce(outer_merge, dfs)
    # Get only LD-independent, common variants
    df = df[(df.wb_maf >= 0.01) & (df.ld_indep == True)]
    # Get rid of unneeded rows
    df = df[
        conserved_columns
        + [
            col
            for col in df.columns
            if (("BETA" in col) or ("OR" in col) or ("P" in col))
        ]
    ]
    df = df.dropna(axis=1, how="all")
    df = df.dropna()
    err_corr = build_err_corr(S, K, pops, phenos, df)
    R_phen = build_R_phen(S, K, pops, phenos, df)
    return err_corr, R_phen


def merge_with_metadata(file_paths, sumstat_files):

    """
    Merges summary statistic files with associated metadata.

    Parameters:
    file_paths: Paths to summary statistic files.
    sumstat_files: Corresponding summary statistic files.

    Returns:
    metadata_sumstat_files: Summary statistic files with additional columns denoting:
        - Gene symbol
        - Most severe consequence
        - MAF in white british population
        - LD independence
    """

    metadata_sumstat_files = []
    exome_metadata = pd.read_csv(
        "/oak/stanford/groups/mrivas/ukbb24983/exome/pgen/spb/data/ukb_exm_spb-consequence_wb_maf_gene_ld_indep.tsv",
        sep="\t",
    )
    cal_metadata = pd.read_csv(
        "/oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb_cal-consequence_wb_maf_gene_ld_indep.tsv",
        sep="\t",
    )
    for file_path, sumstat_file in zip(file_paths, sumstat_files):
        # Merge metadata
        if "exome" in file_path:
            sumstat_file = sumstat_file.merge(exome_metadata)
        elif "cal" in file_path:
            sumstat_file = sumstat_file.merge(cal_metadata)
        metadata_sumstat_files.append(sumstat_file)
    return metadata_sumstat_files


def rename_columns(df, conserved_columns, pop, pheno):

    """ 
    Renames columns such that information on population/study and phenotype is available 
        in the resultant dataframe.
  
    Additionally checks if the header contains "LOG(OR)_SE" instead of "SE".
  
    Parameters: 
    df: Input dataframe (from summary statistics).
    conserved_columns: These columns are not renamed so that we can merge summary 
        statistics across phenotypes and studies.
    pop: The study from which the current summary statistic dataframef comes from.
    pheno: The phenotype from which the current summary statistic dataframe comes from.
  
    Returns: 
    df: A df with adjusted column names, e.g., "OR_white_british_cancer1085".
  
    """

    if "LOG(OR)_SE" in df.columns:
        df.rename(columns={"LOG(OR)_SE": "SE"}, inplace=True)
    columns_to_rename = list(set(df.columns) - set(conserved_columns))
    renamed_columns = [(x + "_" + pop + "_" + pheno) for x in columns_to_rename]
    df.rename(columns=dict(zip(columns_to_rename, renamed_columns)), inplace=True)
    return df


def read_in_summary_stats(pops, phenos, dataset, conserved_columns):

    """ 
    Reads in GBE summary statistics from the Rivas Lab file organization system on 
        Sherlock.
  
    Additionally: adds a variant identifier ("V"), renames columns, and filters on 
        SE (<= 0.5).

    Contains logic for handling the case that a summary statistic file is not found.
  
    Parameters: 
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    dataset: One of ("cal"/"exome") to use for analysis.
    conserved_columns: Columns in every annotated summary statistic file; 
        basis of merges.
  
    Returns: 
    file_paths: List of strings containing file paths.
    sumstat_files: List of dataframes containing summary statistics.
  
    """

    file_paths = []
    sumstat_files = []
    for pop in pops:
        for pheno in phenos:
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


def collect_and_filter(pops, phenos, dataset, conserved_columns, maf_thresh):

    """ 
    Collects the summary statistics of interest and applies filters. 
  
    Reads in summary statistics; 
    Joins on metadata files that contain gene symbol, MAF, and consequence;
    Sets sigma values associated with different types of consequences (PTV/PAV/PCV).
    
    Parameters: 
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    dataset: One of ("cal"/"exome") to use for analysis.
    conserved_columns: Columns in every annotated summary statistic file; basis of merges.
    maf_thresh: Maximum MAF of variants in this run.
  
    Returns: 
    filtered_sumstats_files: Annotated and filtered summary statistic files.
    err_corr: (S*K x S*K) matrix of correlation of errors across studies and phenotypes 
        for null variants. Used to calculate v_beta.
  
    """

    print(Fore.CYAN + "Reading in summary statistics for:" + Style.RESET_ALL)
    print("")
    print(Fore.CYAN + "Populations: " + Style.RESET_ALL + ", ".join(pops))
    print(Fore.CYAN + "Phenotypes: " + Style.RESET_ALL + ", ".join(phenos))
    print(Fore.CYAN + "Dataset: " + Style.RESET_ALL + dataset)
    print("")
    file_paths, sumstat_files = read_in_summary_stats(
        pops, phenos, dataset, conserved_columns
    )
    metadata_sumstat_files = merge_with_metadata(file_paths, sumstat_files)
    err_corr, R_phen = return_err_and_R_phen(
        metadata_sumstat_files, pops, phenos, len(pops), len(phenos)
    )
    filtered_sumstat_files = filter_and_set_sigmas(metadata_sumstat_files, maf_thresh)
    return filtered_sumstat_files, err_corr, R_phen


def return_input_args(args):

    """ 
    Further parses the command-line input.
  
    Makes all lists unique; calculates S and K; and creates lists of appropriate 
        matrices.
  
    Parameters: 
    args: Command-line arguments that have been parsed by the parser.
  
    Returns: 
    S: Number of populations/studies.
    K: Number of GBE phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of GBE phenotypes to use for analysis.
    R_study_list: Unique list of R_study matrices to use for analysis.
    R_study_models: Unique strings ("independent"/"similar") corresponding to each 
        matrix in R_study_list.
    R_var_models: Unique strings ("independent"/"similar") corresponding to R_var 
        matrices to use for analysis.
    agg: Unique list of aggregation units ("gene"/"variant") to use for analysis.
    sigma_m_types: Unique list of sigma_m types ("sigma_m_var"/"sigma_m_1"/"sigma_m_005") 
        to use for analysis.
    variant_filters: Unique list of variant filters ("ptv"/"pav"/"pcv") to use 
        for analysis.
    datasets: Unique list of datasets ("cal"/"exome") to use for analysis.
    maf_threshes: List of maximum MAFs of variants in your runs.
    prior_odds_list: List of prior odds used as assumptions to calculate posterior 
        probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from 
        Bayes Factors.
  
    """

    for arg in vars(args):
        setattr(args, arg, list(set(getattr(args, arg))))
    S = len(args.pops)
    K = len(args.phenos)
    R_study = [
        np.diag(np.ones(S)) if x == "independent" else np.ones((S, S))
        for x in args.R_study_models
    ]
    return (
        S,
        K,
        args.pops,
        args.phenos,
        R_study,
        args.R_study_models,
        args.R_var_models,
        args.agg,
        args.sigma_m_types,
        args.variants,
        args.datasets,
        args.maf_threshes,
        args.prior_odds_list,
        args.p_value_methods,
    )


def print_banner():

    """ 
    Prints ASCII Art Banner + Author Info.
  
    """

    print(Fore.RED + " __  __ ____  ____")
    print("|  \/  |  _ \|  _ \\")
    print("| |\/| | |_) | |_) |")
    print("| |  | |  _ <|  __/ ")
    print("|_|  |_|_| \_\_|  " + Style.RESET_ALL)
    print("")
    print(Fore.GREEN + "Production Author:" + Style.RESET_ALL)
    print("Guhan Ram Venkataraman, B.S.H.")
    print("Ph.D. Candidate | Biomedical Informatics")
    print("")
    print(Fore.GREEN + "Contact:" + Style.RESET_ALL)
    print("Email: guhan@stanford.edu")
    print(
        "URL: https://github.com/rivas-lab/ukbb-tools/blob/master/13_mrp/mrp_production.py"
    )
    print("")
    print(Fore.GREEN + "Methods Developers:" + Style.RESET_ALL)
    print("Manuel A. Rivas, Ph.D.; Matti Pirinen, Ph.D.")
    print("Rivas Lab | Stanford University")
    print("")


def range_limited_float_type(arg):

    """ 
    Type function for argparse - a float within some predefined bounds.
    
    Parameters:
    arg: Putative float.

    """

    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "--maf_thresh must be a valid floating point number."
        )
    if f <= 0 or f >= 1:
        raise argparse.ArgumentTypeError("--maf_thresh must be > 0 and < 1.")
    return f


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
        help="name(s) of populations/studies to be meta-analyzed; at least one required.",
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
        help="name(s) of phenotypes to be studied jointly; at least one required.",
    )
    parser.add_argument(
        "--R_var",
        choices=["independent", "similar"],
        type=str,
        nargs="+",
        default=["independent"],
        dest="R_var_models",
        help="type(s) of model across variants. options: independent, similar (default: independent). can run both.",
    )
    parser.add_argument(
        "--M",
        choices=["variant", "gene"],
        type=str,
        nargs="+",
        default=["gene"],
        dest="agg",
        help="unit(s) of aggregation. options: variant, gene (default: gene). can run both.",
    )
    parser.add_argument(
        "--sigma_m_types",
        choices=["sigma_m_var", "sigma_m_1", "sigma_m_005"],
        type=str,
        nargs="+",
        default=["sigma_m_var"],
        dest="sigma_m_types",
        help="scaling factor(s) for variants. options: var (i.e. 0.2 for ptvs, 0.05 for pavs/pcvs), 1, 0.05 (default: var). can run multiple.",
    )
    parser.add_argument(
        "--variants",
        choices=["pcv", "pav", "ptv"],
        type=str,
        nargs="+",
        default=["ptv"],
        dest="variants",
        help="variant set(s) to consider. options: proximal coding [pcv], protein-altering [pav], protein truncating [ptv] (default: ptv). can run multiple.",
    )
    parser.add_argument(
        "--datasets",
        choices=["cal", "exome"],
        type=str,
        nargs="+",
        default=["cal"],
        dest="datasets",
        help="which UKBB dataset(s) to use. options: cal, exome (default: cal). can run multiple.",
    )
    parser.add_argument(
        "--maf_thresh",
        type=range_limited_float_type,
        nargs="+",
        default=[0.01],
        dest="maf_threshes",
        help="which MAF threshold(s) to use. must be valid floats between 0 and 1 (default: 0.01).",
    )
    parser.add_argument(
        "--prior_odds",
        type=range_limited_float_type,
        nargs="+",
        default=[0.0005],
        dest="prior_odds_list",
        help="which prior odds (can be multiple) to use in calculating posterior probabilities. must be valid floats between 0 and 1 (default: 0.0005, expect 1 in 2000 genes to be a discovery).",
    )
    parser.add_argument(
        "--p_value",
        choices=["farebrother", "davies", "imhof"],
        type=str,
        nargs="+",
        default=[],
        dest="p_value_methods",
        help="which method(s) to use to convert Bayes Factors to p-values. if command line argument is invoked but method is not specified, will throw an error (i.e., specify a method when it is invoked). if not invoked, p-values will not be calculated. options: farebrother, davies, imhof. NOTE: --p_value imports R objects and methods, which slows down MRP. farebrother is fastest and recommended if p-values are a must.",
    )
    return parser


if __name__ == "__main__":

    """ 
    Runs MRP analysis on GBE summary statistics with the parameters specified by the command line.

    """

    with open("../05_gbe/phenotype_info.tsv", "r") as phe_file:
        valid_phenos = [line.split()[0] for line in phe_file][1:]
    parser = initialize_parser(valid_phenos)
    args = parser.parse_args()
    print("")
    print("Valid arguments. Importing required packages...")
    print("")
    import pandas as pd
    from functools import partial, reduce

    pd.options.mode.chained_assignment = None
    import numpy as np
    import numpy.matlib as npm
    from numpy.linalg import LinAlgError
    from scipy.stats.stats import pearsonr
    import subprocess
    import os
    from colorama import Fore, Back, Style

    print_banner()
    S, K, pops, phenos, R_study_list, R_study_models, R_var_models, agg, sigma_m_types, variant_filters, datasets, maf_threshes, prior_odds_list, p_value_methods = return_input_args(
        args
    )
    if p_value_methods:
        print(
            Fore.RED
            + "WARNING: Command line arguments indicate p-value generation. This can cause slowdowns of up to 12x."
        )
        print(
            "Consider using --prior_odds instead to generate posterior probabilities as opposed to p-values."
            + Style.RESET_ALL
        )
        print("")
        import rpy2
        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri

        rpy2.robjects.numpy2ri.activate()
        import rpy2.robjects.packages as rpackages
        from rpy2.robjects.vectors import StrVector
        from rpy2.robjects.vectors import ListVector
        from rpy2.robjects.vectors import FloatVector
        import warnings
        from rpy2.rinterface import RRuntimeWarning

        warnings.filterwarnings("ignore", category=RRuntimeWarning)
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
    for maf_thresh in maf_threshes:
        for dataset in datasets:
            sumstat_files, err_corr, R_phen = collect_and_filter(
                pops, phenos, dataset, conserved_columns, maf_thresh
            )
            loop_through_parameters(
                dataset,
                agg,
                variant_filters,
                S,
                R_study_list,
                R_study_models,
                pops,
                K,
                R_phen,
                phenos,
                R_var_models,
                sigma_m_types,
                sumstat_files,
                err_corr,
                conserved_columns,
                maf_thresh,
                prior_odds_list,
                p_value_methods,
            )