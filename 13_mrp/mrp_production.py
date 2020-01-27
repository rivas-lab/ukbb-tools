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
    # Check that it's not only pos def but also above a certain threshold
    X = np.matrix(X)
    i = 0
    while (np.linalg.cond(X) >= 1/np.finfo(X.dtype).eps):
        X = 0.99 * X + 0.01 * np.diag(np.diag(X))
        i += 1
        if i >= 5:
            return X
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
        X = is_pos_def_and_full_rank(X)
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
    if np.any(np.isnan(A)):
        return [np.nan] * len(methods)
    A_inv = np.linalg.inv(A)
    quad_T = np.asmatrix(beta.T) * np.asmatrix((v_beta_inv - A_inv)) * np.asmatrix(beta)
    B = is_pos_def_and_full_rank(npm.eye(n) - np.asmatrix(A_inv) * np.asmatrix(v_beta))
    if np.any(np.isnan(B)):
        return [np.nan] * len(methods)
    d = np.linalg.eig(B)[0]
    d = [i for i in d if i > 0.01]
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
    posterior_probs: List of posterior probabilities of the event
        given the list of prior odds and the Bayes Factor.
  
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
    v_beta_inv = safe_inv(v_beta, "v_beta", block, agg_type)
    U_inv = safe_inv(U, "U", block, agg_type)
    if v_beta_inv is not np.nan and U_inv is not np.nan:
        A2 = U_inv + v_beta_inv
        b2 = np.asmatrix(v_beta_inv) * np.asmatrix(beta)
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
    phenos: Unique set of phenotypes to use for analysis.
  
    Returns: 
    beta_list: A list of effect sizes (some may be missing) from the subset.
    se_list: A list of standard errors (some may be missing) from the subset.
  
    """

    beta_list = []
    se_list = []
    for pop in pops:
        for pheno in phenos:
            beta_list.extend(list(subset_df["BETA" + "_" + pop + "_" + pheno]))
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
    phenos: Unique set of phenotypes to use for analysis.
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
    print("diag_sigma_m")
    print(diag_sigma_m)
    np.save("HADH_diag_sigma_m", diag_sigma_m)
    R_var = np.diag(np.ones(M)) if R_var_model == "independent" else np.ones((M, M))
    S_var = np.dot(np.dot(diag_sigma_m, R_var), diag_sigma_m)
    if R_var_model == "independent":
        np.save("HADH_Svar_independent.npy", S_var)
    else:
        np.save("HADH_Svar_similar.npy", S_var)
    beta_list, se_list = generate_beta_se(subset_df, pops, phenos)
    beta = np.array(beta_list).reshape(-1, 1)
    se = np.array(se_list)
    omega = np.kron(err_corr, np.diag(np.ones(M)))
    U = np.kron(np.kron(R_study, R_phen), S_var)
    U, omega, beta, se = adjust_for_missingness(U, omega, beta, se, beta_list)
    print("se")
    print(se)
    np.save('HADH_se.npy', se)
    diag_se = np.diag(se)
    v_beta = np.dot(np.dot(diag_se, omega), diag_se)
    mu = np.zeros(beta.shape)
    return U, beta, v_beta, mu


def output_file(bf_dfs, agg_type, pops, phenos, maf_thresh, out_folder):

    """ 
    Outputs a file containing aggregation unit and Bayes Factors. 
    
    Parameters: 
    bf_dfs: List of dataframes containing Bayes Factors from each analysis.
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    maf_thresh: Maximum MAF of variants in this run.
    out_folder: Output folder in which results are stored.

    """

    outer_merge = partial(pd.merge, on=agg_type, how="outer")
    out_df = reduce(outer_merge, bf_dfs)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
        print("")
        print(Fore.RED + "Folder " + out_folder + " created." + Style.RESET_ALL)
        print("")
    out_file = os.path.join(
        out_folder,
        "_".join(pops)
        + "_"
        + "_".join(phenos)
        + "_"
        + agg_type
        + "_"
        + str(maf_thresh)
        + ".tsv",
    )
    out_df = out_df.sort_values(
        by=out_df.columns[1],
        ascending=False,
    )
    out_df.to_csv(out_file, sep="\t", index=False)
    print("")
    print(Fore.RED + "Results written to " + out_file + "." + Style.RESET_ALL)
    print("")


def get_output_file_columns(
    agg_type,
    R_study_model,
    R_var_model,
    sigma_m_type,
    analysis,
    prior_odds_list,
    p_value_methods,
):

    """
    Sets up the columns that must be enumerated in the output dataframe from MRP.

    Parameters:
    agg_type:
    R_study_model: String ("independent"/"similar") corresponding to R_study.
    R_var_model: String ("independent"/"similar") corresponding to R_var matrices to 
        use for analysis.
    sigma_m_type: One of "sigma_m_var"/"sigma_m_1"/"sigma_m_005". Dictates variant 
        scaling factor by functional annotation.
    analysis: One of "ptv"/"pav"/"pcv". Dictates which variants are included.
    prior_odds_list: List of prior odds used as assumptions to calculate posterior 
        probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from Bayes 
        Factors.

    Returns:
    bf_df_columns: Columns needed for the output file.
    fb: Farebrother R method (rpy2 object), or None if --p_value is not invoked.
    dm: Davies R method (rpy2 object), or None if --p_value is not invoked.
    im: Imhof R method (rpy2 object), or None if --p_value is not invoked.
    
    """

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
    return bf_df_columns, fb, dm, im


def run_mrp(
    df,
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
    agg_type,
    prior_odds_list,
    p_value_methods,
):

    """ 
    Runs MRP with the given parameters.
  
    Parameters: 
    df: Merged dataframe containing all relevant summary statistics.
    S: Number of populations/studies.
    K: Number of phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
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
    agg_type: One of "gene"/"variant". Dictates block of aggregation.
    prior_odds_list: List of prior odds used as assumptions to calculate posterior 
        probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from Bayes 
        Factors.

    Returns: 
    bf_df: Dataframe with two columns: agg_type and log_10 Bayes Factor. 
  
    """
    m_dict = (
        df.groupby("gene_symbol").size()
        if agg_type == "gene"
        else df.groupby("V").size()
    )
    bf_df_columns, fb, dm, im = get_output_file_columns(
        agg_type,
        R_study_model,
        R_var_model,
        sigma_m_type,
        analysis,
        prior_odds_list,
        p_value_methods,
    )
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
        print("U")
        print(U)
        np.save("HADH_U.npy", U)
        print("beta")
        print(beta)
        np.save("HADH_beta.npy", beta)
        print("v_beta")
        print(v_beta)
        np.save("HADH_v_beta.npy", v_beta)
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


def filter_category(df, variant_filter):

    """ 
    Filters a set of dataframes that have been read in based on functional consequence.
  
    Dependent on the variant filter that is dictated by the analysis.
  
    Parameters: 
    df: Merged dataframe containing all summary statistics.
    variant_filter: The variant filter dictated by the analysis ("ptv"/"pav"/"pcv").
  
    Returns: 
    df: Merged dataframe containing all relevant summary statistics; 
        filters out variants excluded from analysis.
  
    """

    if variant_filter == "ptv":
        df = df[df.category == "ptv"]
    elif variant_filter == "pav":
        df = df[(df.category == "ptv") | (df.category == "pav")]
    return df


def loop_through_parameters(
    df,
    maf_threshes,
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
    err_corr,
    prior_odds_list,
    p_value_methods,
    out_folder,
):

    """ 
    Loops through parameters specified through command line (or defaults). 

    Parameters: 
    df: Merged dataframe containing all summary statistics.
    maf_threshes: List of maximum MAFs of variants in your runs.
    agg: Unique list of aggregation units ("gene"/"variant") to use for analysis.
    variant_filters: Unique list of variant filters ("ptv"/"pav"/"pcv") to use 
        for analysis.
    S: Number of populations/studies.
    R_study_list: Unique list of R_study matrices to use for analysis.
    R_study_models: Unique strings ("independent"/"similar") corresponding to each 
        matrix in R_study_list.
    pops: Unique set of populations (studies) to use for analysis.
    K: Number of phenotypes.
    R_phen: R_phen matrix to use for analysis (empirically calculated).
    phenos: Unique set of phenotypes to use for analysis.
    R_var_models: Unique strings ("independent"/"similar") corresponding to R_var 
        matrices to use for analysis.
    sigma_m_types: Unique list of sigma_m types ("sigma_m_var"/"sigma_m_1"/"sigma_m_005")
        to use for analysis.
    err_corr: Matrix of correlation of errors across studies and phenotypes.
    prior_odds_list: List of prior odds used as assumptions to calculate posterior 
        probabilities of Bayes Factors.
    p_value_methods: List of p-value methods used to calculate p-values from Bayes 
        Factors.
    out_folder: Folder where output will be placed.
  
    """

    if (S == 1) and (len(R_study_models) > 1):
        print(
            Fore.YELLOW
            + "Since we are not meta-analyzing, R_study is just [1]."
            + Style.RESET_ALL
        )
        print("")
        R_study_models = ["independent"]
        R_study_list = [R_study_list[0]]
    for maf_thresh in maf_threshes:
        print(
            Fore.YELLOW
            + "Running MRP across parameters for maf_thresh "
            + str(maf_thresh)
            + "..."
            + Style.RESET_ALL
        )
        maf_df = df[(df.maf <= maf_thresh) & (df.maf > 0)]
        for agg_type in agg:
            bf_dfs = []
            # If not aggregating, then R_var choice does not affect BF
            if (agg_type == "variant") and (len(R_var_models) > 1):
                print(
                    Fore.YELLOW
                    + "Since we are not aggregating, R_var is just [1]."
                    + Style.RESET_ALL
                )
                R_var_models = ["independent"]
            for analysis in variant_filters:
                analysis_df = filter_category(maf_df, analysis)
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
                                analysis_df,
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
                                agg_type,
                                prior_odds_list,
                                p_value_methods,
                            )
                            bf_dfs.append(bf_df)
            output_file(bf_dfs, agg_type, pops, phenos, maf_thresh, out_folder)


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
        return [], []
    if mode == "null":
        df = df[
            (df["P_" + pop1 + "_" + pheno1].astype(float) >= 1e-2)
            & (df["P_" + pop2 + "_" + pheno2].astype(float) >= 1e-2)
        ]
    elif mode == "sig":
        df = df[
            ((df["P_" + pop1 + "_" + pheno1].astype(float) <= 1e-5)
            | (df["P_" + pop2 + "_" + pheno2].astype(float) <= 1e-5))
        ]
    beta1 = list(df["BETA_" + pop1 + "_" + pheno1])
    beta2 = list(df["BETA_" + pop2 + "_" + pheno2])
    return beta1, beta2


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
            phen_beta1, phen_beta2 = get_betas(df, pop1, pheno1, pop2, pheno2, "sig")
            return (
                pearsonr(phen_beta1, phen_beta2)[0]
                if phen_beta1 is not None
                else np.nan
            )
        else:
            return np.nan


def build_phen_corr(S, K, pops, phenos, df, pop_pheno_tuples):

    """
    Builds out a matrix of correlations between all phenotypes and studies using:
        - significant (P < 1e-5)
        - common (MAF >= 0.01)
        - LD-independent
    SNPs.

    Parameters:
    S: Number of populations/studies.
    K: Number of phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.
    pop_pheno_tuples: Indicate which populations/phenotypes to use to build R_phen.

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
                    phen_corr[a, b] = calculate_phen(
                        a, b, pop1, pheno1, pop2, pheno2, df, pop_pheno_tuples
                    )
    return phen_corr


def filter_for_phen_corr(df, map_file):

    """
    Filters the initial dataframe for the criteria used to build R_phen.

    Parameters:
    df: Merged dataframe containing all summary statistics.
    map_file: Dataframe indicating which summary statistics to use to build R_phen.

    Returns:
    df: Filtered dataframe that contains significant, common, LD-independent variants.

    """

    files_to_use = map_file[map_file["R_phen"] == True]
    if len(files_to_use) == 0:
        return [], []
    pop_pheno_tuples = zip(list(files_to_use["study"]), list(files_to_use["pheno"]))
    cols_to_keep = ["V", "maf", "ld_indep"]
    for col_type in "BETA_", "P_":
        cols_to_keep.extend(
            [col_type + pop + "_" + pheno for pop, pheno in pop_pheno_tuples]
        )
    df = df[cols_to_keep]
    # Get only LD-independent, common variants
    df = df[(df.maf >= 0.01) & (df.ld_indep == True)]
    df = df.dropna(axis=1, how="all")
    df = df.dropna()
    return df, pop_pheno_tuples


def build_R_phen(S, K, pops, phenos, df, map_file):

    """
    Builds R_phen using phen_corr (calculated using the method directly above this).

    Parameters:
    S: Number of populations/studies.
    K: Number of phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.
    map_file: Input file containing summary statistic paths + pop and pheno data.

    Returns:
    R_phen: Empirical estimates of genetic correlation across phenotypes.

    """

    if K == 1:
        return np.ones((K, K))
    df, pop_pheno_tuples = filter_for_phen_corr(df, map_file)
    if len(df) == 0:
        print("")
        print(Fore.RED + "WARNING: No files specified for R_phen generation.")
        print("Assuming independent effects." + Style.RESET_ALL)
        return np.diag(np.ones(K))
    phen_corr = build_phen_corr(S, K, pops, phenos, df, pop_pheno_tuples)
    R_phen = np.zeros((K, K))
    for k1, pheno1 in zip(range(K), phenos):
        for k2, pheno2 in zip(range(K), phenos):
            if k1 == k2:
                R_phen[k1, k2] = 1
            elif k1 > k2:
                R_phen[k1, k2] = R_phen[k2, k1]
            else:
                phenos_to_remove = list(set(range(K)) - set([k1, k2]))
                indices_to_remove = []
                for pheno_to_remove in phenos_to_remove:
                    indices_to_remove.extend(
                        [pheno_to_remove + K * x for x in range(S)]
                    )
                pairwise_corrs = delete_rows_and_columns(phen_corr, indices_to_remove)
                R_phen[k1, k2] = np.nanmedian(pairwise_corrs)
    R_phen = np.nan_to_num(R_phen)
    return R_phen


def calculate_err(a, b, pop1, pheno1, pop2, pheno2, err_corr, err_df):

    """
    Calculates a single entry in the err_corr matrix.
    
    Parameters:
    a, b: Positional parameters within the err_corr matrix.
    pop1: Name of first population.
    pheno1: Name of first phenotype.
    pop2: Name of second population.
    pheno2: Name of second phenotype.
    err_corr: The err_corr matrix thus far.
    err_df: Dataframe containing null, common, LD-independent variants.

    Returns:
    err_corr[a, b]: One entry in the err_corr matrix.

    """

    # If in lower triangle, do not compute; symmetric matrix
    if a > b:
        return err_corr[b, a]
    elif a == b:
        return 1
    else:
        err_df = err_df.dropna()
        err_beta1, err_beta2 = get_betas(err_df, pop1, pheno1, pop2, pheno2, "null")
        return pearsonr(err_beta1, err_beta2)[0] if err_beta1 else 0


def filter_for_err_corr(df):

    """
    Filters the initial dataframe for the criteria used to build err_corr.

    Parameters:
    df: Merged dataframe containing all summary statistics.

    Returns:
    df: Filtered dataframe that contains null, common, LD-independent variants.

    """

    print("")
    print(
        Fore.MAGENTA
        + "Building R_phen and matrix of correlations of errors..."
        + Style.RESET_ALL
    )
    print("")
    # Get only LD-independent, common variants
    df = df[(df.maf >= 0.01) & (df.ld_indep == True)]
    df = df.dropna(axis=1, how="all")
    df = df.dropna()
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
    if len(df) != 0:
        df = df[df.most_severe_consequence.isin(null_variants)]
    return df


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
    K: Number of phenotypes.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.

    Returns:
    err_corr: (S*K x S*K) matrix of correlation of errors across studies and phenotypes 
        for null variants. Used to calculate v_beta.

    """
    err_df = filter_for_err_corr(df)
    if len(err_df) == 0:
        print(Fore.RED + "WARNING: Correlation of errors is noisy.")
        print("Assuming independent effects." + Style.RESET_ALL)
        print("")
        return np.diag(np.ones(S * K))
    err_corr = np.zeros((S * K, S * K))
    for i, pop1 in enumerate(pops):
        for j, pheno1 in enumerate(phenos):
            for x, pop2 in enumerate(pops):
                for y, pheno2 in enumerate(phenos):
                    # Location in matrix
                    a, b = K * i + j, K * x + y
                    err_corr[a, b] = calculate_err(
                        a, b, pop1, pheno1, pop2, pheno2, err_corr, err_df
                    )
    err_corr = np.nan_to_num(err_corr)
    return err_corr


def return_err_and_R_phen(df, pops, phenos, S, K, map_file):

    """ 
    Builds a matrix of correlations of errors across studies and phenotypes,
        and correlations of phenotypes.
  
    Parameters: 
    df: Dataframe that containa summary statistics.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    S: Number of populations/studies.
    K: Number of phenotypes.
    map_file: Input file containing summary statistic paths + pop and pheno data.

    Returns:
    err_corr: (S*K x S*K) matrix of correlation of errors across studies and phenotypes
        for null variants. Used to calculate v_beta.
    R_phen: Empirical estimates of genetic correlation across phenotypes.
  
    """

    # Sample common variants, stuff in filter + synonymous
    err_corr = build_err_corr(S, K, pops, phenos, df)
    # Faster calculations, better accounts for uncertainty in estimates
    err_corr[abs(err_corr) < 0.01] = 0
    R_phen = build_R_phen(S, K, pops, phenos, df, map_file)
    R_phen[abs(R_phen) < 0.01] = 0
    # Get rid of any values above 0.95
    while (np.max(R_phen - np.eye(len(R_phen))) > 0.9):
        R_phen = 0.9 * R_phen + 0.1 * np.diag(np.diag(R_phen))
    return err_corr, R_phen


def rename_columns(df, pop, pheno):

    """ 
    Renames columns such that information on population/study and phenotype is available 
        in the resultant dataframe.
  
    Additionally checks if the header contains "LOG(OR)_SE" instead of "SE".
  
    Parameters: 
    df: Input dataframe (from summary statistics).
    pop: The study from which the current summary statistic dataframef comes from.
    pheno: The phenotype from which the current summary statistic dataframe comes from.
  
    Returns: 
    df: A df with adjusted column names, e.g., "OR_white_british_cancer1085".
  
    """

    if "LOG(OR)_SE" in df.columns:
        df.rename(columns={"LOG(OR)_SE": "SE"}, inplace=True)
    columns_to_rename = ["BETA", "SE", "P"]
    renamed_columns = [(x + "_" + pop + "_" + pheno) for x in columns_to_rename]
    df.rename(columns=dict(zip(columns_to_rename, renamed_columns)), inplace=True)
    return df


def check_map_file(map_file):

    """
    Checks --file for malformed input.

    Parameters:
    map_file: Input file containing summary statistic paths + pop and pheno data.
    
    Returns:
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    S: Number of populations/studies.
    K: Number of phenotypes.

    """
    if map_file.isnull().values.sum() > 0:
        raise ValueError("NaNs in map file.")
    file_paths = np.unique(list(map_file["path"]))
    if len(file_paths) < len(map_file):
        raise ValueError("File specified in map file contains duplicate path entries.")
    for file_path in file_paths:
        if not os.path.exists(file_path):
            raise IOError("File " + file_path + ", listed in map file does not exist.")
    if (map_file.groupby(["study", "pheno"]).size() > 1).sum() > 0:
        raise ValueError(
            "Multiple summary statistic files specified for a (study, phenotype) tuple."
        )
    pops = sorted(list(set(map_file["study"])))
    phenos = sorted(list(set(map_file["pheno"])))
    S = len(pops)
    K = len(phenos)
    return pops, phenos, S, K


def merge_dfs(sumstat_files, metadata_path):

    """
    Performs an outer merge on all of the files that have been read in;
    Annotates with metadata and sigma values.

    Parameters:
    sumstat_files: List of dataframes that contain summary statistics.
    metadata_path: Path to metadata file containing MAF, Gene symbol, etc.

    Returns:
    df: Dataframe that is ready for err_corr/R_phen generation and for running MRP.

    """

    print("")
    print(Fore.CYAN + "Merging summary statistics with metadata...")
    conserved_columns = ["V", "#CHROM", "POS", "REF", "ALT", "A1"]
    outer_merge = partial(pd.merge, on=conserved_columns, how="outer")
    df = reduce(outer_merge, sumstat_files)
    metadata = pd.read_csv(metadata_path, sep="\t")
    df = df.merge(metadata)
    df = set_sigmas(df)
    return df


def read_in_summary_stat(subset_df, pop, pheno):

    """
    Reads in one summary statistics file.
  
    Additionally: adds a variant identifier ("V"), renames columns, and filters on 
        SE (<= 0.5).

    Parameters: 
    subset_df: Subset of the map file where study == pop and phenotype == pheno.
    pop: Population of interest.
    pheno: Phenotype of interest.
  
    Returns: 
    df: Dataframe with renamed columns, ready for merge.

    """

    file_path = list(subset_df["path"])[0]
    print(file_path)
    df = pd.read_csv(
        file_path,
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
    if "OR" in df.columns:
        df["BETA"] = np.log(df["OR"].astype("float64"))
    # Filter for SE as you read it in
    df = rename_columns(df, pop, pheno)
    df = df[df["SE" + "_" + pop + "_" + pheno].notnull()]
    df = df[df["SE" + "_" + pop + "_" + pheno].astype(float) <= 0.2]
    # Filter out HLA region
    df = df[~((df["#CHROM"] == 6) & (df["POS"].between(25477797, 36448354)))]
    return df


def read_in_summary_stats(map_file, metadata_path, exclude_path):

    """ 
    Reads in summary statistics.
  
    Additionally: adds a variant identifier ("V"), renames columns, and filters on 
        SE (<= 0.5).

    Contains logic for handling the case that a summary statistic file is not found.
  
    Parameters: 
    map_file: Input file containing summary statistic paths + pop and pheno data.
    metadata_path: Path to metadata file containing MAF, Gene symbol, etc.
  
    Returns: 
    df: Merged summary statistics.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    S: Number of populations/studies.
    K: Number of phenotypes. 

    """

    pops, phenos, S, K = check_map_file(map_file)
    print(Fore.CYAN + "Map file passes initial checks.")
    print("")
    print("Reading in summary statistics for:")
    print("")
    print("Populations: " + Style.RESET_ALL + ", ".join(pops))
    print(Fore.CYAN + "Phenotypes: " + Style.RESET_ALL + ", ".join(phenos))
    print("")
    sumstat_files = []
    for pop in pops:
        for pheno in phenos:
            subset_df = map_file[(map_file.study == pop) & (map_file.pheno == pheno)]
            if len(subset_df) == 1:
                df = read_in_summary_stat(subset_df, pop, pheno)
                sumstat_files.append(df)
            else:
                print(
                    Fore.RED
                    + "WARNING: A summary statistic file cannot be found for "
                    + "population: {}; phenotype: {}.".format(pop, pheno)
                    + Style.RESET_ALL
                )
    try:
        if exclude_path:
            exclude_path = exclude_path[0]
            variants_to_exclude = [line.rstrip('\n') for line in open(exclude_path)]
    except:
        raise IOError("Could not open exclusions file (--exclude).")
    df = merge_dfs(sumstat_files, metadata_path)
    if exclude_path:
        df = df[~df["V"].isin(variants_to_exclude)]
    return df, pops, phenos, S, K


def print_banner():

    """ 
    Prints ASCII Art Banner + Author Info.
  
    """

    print("")
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


def return_input_args(args):

    """ 
    Further parses the command-line input.
  
    Makes all lists unique; calculates S and K; and creates lists of appropriate 
        matrices.
  
    Parameters: 
    args: Command-line arguments that have been parsed by the parser.
  
    Returns: 
    df: Merged summary statistics.
    pops: Unique set of populations (studies) to use for analysis.
    phenos: Unique set of phenotypes to use for analysis.
    S: Number of populations/studies.
    K: Number of phenotypes. 
    R_study_list: Unique list of R_study matrices to use for analysis.
  
    """

    try:
        map_file = pd.read_csv(args.map_file, sep="\t")
    except:
        raise IOError("File specified in --file does not exist.")
    df, pops, phenos, S, K = read_in_summary_stats(map_file, args.metadata_path, args.exclude)
    for arg in vars(args):
        setattr(args, arg, sorted(list(set(getattr(args, arg)))))
    R_study = [
        np.diag(np.ones(S)) if x == "independent" else np.ones((S, S))
        for x in args.R_study_models
    ]
    return (df, map_file, S, K, pops, phenos, R_study)


def range_limited_float_type(arg):

    """ 
    Type function for argparse - a float within some predefined bounds.
    
    Parameters:
    arg: Putative float.

    """

    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a valid floating point number.")
    if f <= 0 or f >= 1:
        raise argparse.ArgumentTypeError("must be > 0 and < 1.")
    return f


def initialize_parser(valid_phenos):

    """
    Parses inputs using argparse. 
    
    Parameters: 
    valid_phenos: List of valid phenotypes.

    """

    parser = argparse.ArgumentParser(
        description="MRP takes in several variables that affect how it runs.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--file",
        type=str,
        required=True,
        dest="map_file",
        help="""path to tab-separated file containing list of: 
         summary statistic file paths,
         corresponding studies,
         phenotypes, and
         whether or not to use the file in R_phen generation.
       
         format:
         
         path        study        pheno        R_phen
         /path/to/file1   study1    pheno1     TRUE
         /path/to/file2   study2    pheno1     FALSE
         """,
    )
    parser.add_argument(
        "--metadata_path",
        type=str,
        required=True,
        dest="metadata_path",
        help="""path to tab-separated file containing:
         variants,
         gene symbols,
         consequences,
         MAFs,
         and LD independence info.
       
         format:
         
         V       gene_symbol     most_severe_consequence maf  ld_indep
         1:69081:G:C     OR4F5   5_prime_UTR_variant     0.000189471     False
        """,
    )
    parser.add_argument(
        "--R_study",
        choices=["independent", "similar"],
        type=str,
        nargs="+",
        default=["similar"],
        dest="R_study_models",
        help="""type of model across studies. 
         options: independent, similar (default: similar). can run both.""",
    )
    parser.add_argument(
        "--R_var",
        choices=["independent", "similar"],
        type=str,
        nargs="+",
        default=["independent"],
        dest="R_var_models",
        help="""type(s) of model across variants. 
         options: independent, similar (default: independent). can run both.""",
    )
    parser.add_argument(
        "--M",
        choices=["variant", "gene"],
        type=str,
        nargs="+",
        default=["gene"],
        dest="agg",
        help="""unit(s) of aggregation. 
         options: variant, gene (default: gene). can run both.""",
    )
    parser.add_argument(
        "--sigma_m_types",
        choices=["sigma_m_var", "sigma_m_1", "sigma_m_005"],
        type=str,
        nargs="+",
        default=["sigma_m_var"],
        dest="sigma_m_types",
        help="""scaling factor(s) for variants.
         options: var (i.e. 0.2 for ptvs, 0.05 for pavs/pcvs), 
         1, 0.05 (default: var). can run multiple.""",
    )
    parser.add_argument(
        "--variants",
        choices=["pcv", "pav", "ptv"],
        type=str,
        nargs="+",
        default=["ptv"],
        dest="variant_filters",
        help="""variant set(s) to consider. 
         options: proximal coding [pcv], 
                  protein-altering [pav], 
                  protein truncating [ptv] 
                  (default: ptv). can run multiple.""",
    )
    parser.add_argument(
        "--maf_thresh",
        type=range_limited_float_type,
        nargs="+",
        default=[0.01],
        dest="maf_threshes",
        help="""which MAF threshold(s) to use. must be valid floats between 0 and 1 
         (default: 0.01).""",
    )
    parser.add_argument(
        "--prior_odds",
        type=range_limited_float_type,
        nargs="+",
        default=[0.0005],
        dest="prior_odds_list",
        help="""which prior odds (can be multiple) to use in calculating posterior 
         probabilities. must be valid floats between 0 and 1 (default: 0.0005, expect 
         1 in 2000 genes to be a discovery).""",
    )
    parser.add_argument(
        "--p_value",
        choices=["farebrother", "davies", "imhof"],
        type=str,
        nargs="+",
        default=[],
        dest="p_value_methods",
        help="""which method(s) to use to convert Bayes Factors to p-values. if command 
         line argument is invoked but method is not specified, will throw an error 
         (i.e., specify a method when it is invoked). if not invoked, p-values will not 
         be calculated. options: farebrother, davies, imhof. NOTE: --p_value imports R 
         objects and methods, which slows down MRP. farebrother is fastest and 
         recommended if p-values are a must.""",
    )
    parser.add_argument(
        "--exclude",
        type=str,
        nargs=1,
        default=[],
        dest="exclude",
        help="""path to file containing list of variants to exclude from analysis.

         format of file:

         1:69081:G:C
         1:70001:G:A
        """,
    )
    parser.add_argument(
        "--out_folder",
        type=str,
        nargs=1,
        default=[],
        dest="out_folder",
        help="""folder to which output(s) will be written (default: current folder).
         if folder does not exist, it will be created.""",
    )
    return parser


if __name__ == "__main__":

    """ 
    Runs MRP analysis on summary statistics with the parameters specified 
        by the command line.

    """

    with open("/oak/stanford/groups/mrivas/users/guhan/repos/ukbb-tools/05_gbe/phenotype_info.tsv", "r") as phe_file:
        valid_phenos = [line.split()[0] for line in phe_file][1:]
    import os

    parser = initialize_parser(valid_phenos)
    args = parser.parse_args()
    print("")
    print("Valid command line arguments. Importing required packages...")
    print("")
    import pandas as pd
    from functools import partial, reduce

    pd.options.mode.chained_assignment = None
    import numpy as np
    import numpy.matlib as npm
    from numpy.linalg import LinAlgError
    from scipy.stats.stats import pearsonr
    import subprocess
    from colorama import Fore, Back, Style

    df, map_file, S, K, pops, phenos, R_study_list = return_input_args(args)
    out_folder = args.out_folder[0] if args.out_folder else os.getcwd()
    print_banner()
    if args.p_value_methods:
        print("")
        print(
            Fore.RED
            + "WARNING: Command line arguments indicate p-value generation. "
            + "This can cause slowdowns of up to 12x."
        )
        print(
            "Consider using --prior_odds instead to generate posterior probabilities"
            + " as opposed to p-values."
            + Style.RESET_ALL
        )
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
    err_corr, R_phen = return_err_and_R_phen(
        df, pops, phenos, len(pops), len(phenos), map_file
    )
    df = df[(df["gene_symbol"] == "HADH")]# | (df["gene_symbol"] == "ANGPTL7")]
    loop_through_parameters(
        df,
        args.maf_threshes,
        args.agg,
        args.variant_filters,
        S,
        R_study_list,
        args.R_study_models,
        pops,
        K,
        R_phen,
        phenos,
        args.R_var_models,
        args.sigma_m_types,
        err_corr,
        args.prior_odds_list,
        args.p_value_methods,
        out_folder,
    )
