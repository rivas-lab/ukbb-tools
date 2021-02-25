# coding: utf-8
from __future__ import print_function
from __future__ import division
import argparse

# Written by Manuel A. Rivas
# Updated 02.11.2020, Guhan R. Venkataraman


def is_pos_def(X):

    """ 
    Ensures a matrix is positive definite.
  
    Keep diagonals and multiples every other cell by .99.
  
    Parameters: 
    X: Matrix to verify.
  
    Returns: 
    Boolean: Indicator of whether or not X is positive definite.
  
    """

    X = np.matrix(X)
    if np.all(np.linalg.eigvals(X) > 0):
        return True
    else:
        return False


def initialize_MCMC(
    niter, R_phen, R_phen_inv, err_corr, betas, ses, C, R_phen_use, gene_vec, annot_vec
):

    """ 
    Initializes the parameters of the MCMC run.
  
    Parameters: 
    niter: Number of iterations of the run.
    R_phen: K*K matrix of estimated genotype correlations (significant vars).
    R_phen_inv: Inverse of R_phen.
    err_corr: Correlation of errors (common vars).
    betas: An M*K matrix of betas, padded with 0s if missing summary statistics.
    ses: An M*K matrix of standard errors, padded with 0s if missing summary statistics.
    C: Number of hypothesized clusters.
    R_phen_use: Whether or not to use R_phen to initialize parameters.
    gene_vec: Vector of length M of gene symbols for the M variants.
    annot_vec: Vector of length M of functional annotations for the M variants.
  
    Returns: 
    betas: Numpy matrix version of betas.
    ses: Numpy matrix version of standard errors.
    err_corr: Numpy matrix version of correlation of errors.
    K: Number of phenotypes.
    M: Number of variants.
    gene_len: Number of unique genes among M variants.
    annot_len: Number of unique annotations among M variants
    gene_map: List of unique genes.
    annot_map: List of unique annotations.
    alpha: Inverse-gamma prior for pcj.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    bc: Mean effect size per cluster.
    scales: Scale parameters for annotations across clusters.
    delta_m: Indices of cluster memberships for each variant m in gene j.
    maxloglkiter: Max log likelihood of the iteration.
    Theta_0: Prior estimate of genetic correlation across traits; if R_phen_use is
        true, Theta_0 = R_phen. Else, it is the identity matrix.
    Theta_0_inv: Inverse of Theta_0.
    [accept/reject]_mh[1,2,3]: Trackers for MH steps. Step 1 corresponds to updating
        pc, step 2 corresponds to scales, and step 3 to alpha.
    [accept/reject]_mh[1,2,3]_postburnin: Trackers for MH steps post-burn iterations.
  
    """

    print("Running MCMC algorithm with " + str(C) + " clusters...")
    # Convert all parameters to matrices
    betas, ses, err_corr = np.matrix(betas), np.matrix(ses), np.matrix(err_corr)
    # C is the number of clusters, where cluster 1 is the null model cluster
    # Initialize max log likelihood vector - ???????? WHY + 2???
    maxloglkiter = np.zeros((niter + 2, 1))
    # M = number of variants, K = number of phenotypes
    M, K = betas.shape[0], betas.shape[1]
    # Initialize
    # Sigma0 for alternative clusters
    if R_phen_use:
        if is_pos_def(R_phen):
            Theta_0 = R_phen
            Theta_0_inv = R_phen_inv
        else:
            Theta_0 = covariance.shrunk_covariance(R_phen)
            Theta_0_inv = np.linalg.inv(Theta_0)
    else:
        Theta_0 = np.eye(R_phen.shape[0])
        Theta_0_inv = np.linalg.inv(Theta_0)
    # Initialize scale matrix for types of annotations
    gene_set, annot_set = set(gene_vec), set(annot_vec)
    gene_map, annot_map = list(gene_set), list(annot_set)
    gene_len, annot_len = len(gene_set), len(annot_set)
    # For Metropolois Hastings sub-step : keep track of acceptance rate
    accept_mh1, reject_mh1, accept_mh1_postburnin, reject_mh1_postburnin = 0, 0, 0, 0
    accept_mh3, reject_mh3, accept_mh3_postburnin, reject_mh3_postburnin = 0, 0, 0, 0
    accept_mh2, reject_mh2, accept_mh2_postburnin, reject_mh2_postburnin = (
        [0] * annot_len,
        [0] * annot_len,
        [0] * annot_len,
        [0] * annot_len,
    )
    # initialize \alpha : sharing of clusters across genes
    alpha = np.zeros((niter + 2, 1))
    alpha[0, :] = invgamma.rvs(1, 0, 1, size=1)
    # initialize pc (proportions across all variants)
    # store the probabilities (proportions) of cluster memberships
    pc = np.zeros((niter + 2, 1, C))
    pc[0, 0, :] = np.random.dirichlet([1] * C)
    # initialize pcj (proportions for each gene j)
    # store the probabilities (proportions) of cluster memberships for each gene
    pcj = np.zeros((niter + 2, gene_len, C))
    for gene_idx in range(0, gene_len):
        pcj[0, gene_idx, :] = np.random.dirichlet(alpha[0, 0] * pc[0, 0, :])
    # store the mean trait value across the clusters for individuals that are members
    bc = np.zeros((niter + 2, C, K))
    bc[0, 0, :] = np.array([0] * K)
    for c in range(1, C):
        bc[0, c, :] = np.random.multivariate_normal(np.array([0] * K).T, Theta_0)
    scales = np.zeros((niter + 2, annot_len))
    for scaleidx in range(0, annot_len):
        scales[0, scaleidx] = np.power(0.2, 2)
    # initialize variant membership across clusters for each iteration
    delta_m = np.zeros((niter + 2, M))
    return (
        betas,
        ses,
        err_corr,
        K,
        M,
        gene_len,
        annot_len,
        gene_map,
        annot_map,
        alpha,
        pc,
        pcj,
        bc,
        scales,
        delta_m,
        maxloglkiter,
        Theta_0,
        Theta_0_inv,
        accept_mh1,
        accept_mh1_postburnin,
        reject_mh1,
        reject_mh1_postburnin,
        accept_mh2,
        accept_mh2_postburnin,
        reject_mh2,
        reject_mh2_postburnin,
        accept_mh3,
        accept_mh3_postburnin,
        reject_mh3,
        reject_mh3_postburnin,
    )


def return_norm_const(mult, proposal, epsilon):

    """
    Returns the log of the normalizing constant D(z) given the proposal.

                                  C       
                                 ___      
                                 \        
                           Γ  ⋅  /    z_c 
                                 ‾‾‾      
                                c = 1     
             D(z) = ────────────────────────────
                             C            
                           ━┳┳━           
                            ┃┃   Γ ⋅ (z_c)
                           c = 1          

    Parameters:
    mult: Multiplicative factor in conditional. E.g. gamma / alpha.
    proposal: z_c.
    epsilon: Tolerance.

    Returns:
    norm_const: Normalization constant as described above.

    """

    norm_const = math.lgamma(np.sum([mult * i for i in proposal])) - np.sum(
        [math.lgamma(max(mult * i, epsilon)) for i in proposal]
    )
    return norm_const


def return_product_density(mult, proposal, previous, C):

    """
    The density at point x given the normalization constant as defined in
    return_norm_const is:

                                   C               
                                 ____     (z_c - 1)
             p_dir(x|z) = D(z) .  ||   x_c         
                                 c = 1
 
    This function returns the product part of the density (after D(z)).

    Parameters:
    mult: Multiplicative factor in conditional. E.g. gamma / alpha.
    proposal: z_c.
    previous: x_c.
    C: Number of hypothesized clusters.

    Returns:
    density_prod: Product part of the density.

    """

    density_prod = np.sum(
        [(mult * proposal[i] - 1) * np.log(previous[i]) for i in range(0, C)]
    )
    return density_prod


def calculate_l_pdir_num(
    gamma, C, pc_proposal, epsilon, iteration, gene_len, alpha, pc, pcj
):

    """
    Calculates the log of the numerator of the lambda value as described in
    calculate_l_pdir.

    Parameters:
    gamma: Multiplicative part of the proposal value.
    C: Number of hypothesized clusters.
    pc_proposal: Sampled proposal value from the Dirichlet distribution using alpha and
        previous iteration as the parameters.
    epsilon: Tolerance.
    iteration: Iteration number.
    gene_len: Number of unique genes among M variants.
    alpha: Prior for pcj.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.

    Returns:
    l_pdir_num: The log of the numerator of the lambda value as described in
    calculate_l_pdir.

    """

    lhs_num_norm_const = return_norm_const(gamma, pc_proposal, epsilon)
    num_density_prod = return_product_density(
        gamma, pc_proposal, pc[iteration - 1, 0, :], C
    )
    l_pdir_prop = lhs_num_norm_const + num_density_prod
    l_pdir_prop_gene = 0
    rhs_num_norm_const = return_norm_const(
        alpha[iteration - 1, 0], pc_proposal, epsilon
    )
    for gene_idx in range(0, gene_len):
        num_gene_density_prod = return_product_density(
            alpha[iteration - 1, 0], pc_proposal, pcj[iteration - 1, gene_idx, :], C
        )
        l_pdir_prop_gene += num_gene_density_prod + rhs_num_norm_const
    l_pdir_num = l_pdir_prop + l_pdir_prop_gene
    return l_pdir_num


def calculate_l_pdir_den(
    gamma, C, pc_proposal, epsilon, iteration, gene_len, alpha, pc, pcj
):

    """
    Calculates the log of the denominator of the lambda value as described in
    calculate_l_pdir.

    Parameters:
    gamma: Multiplicative part of the proposal value.
    C: Number of hypothesized clusters.
    pc_proposal: Sampled proposal value from the Dirichlet distribution using alpha and
        previous iteration as the parameters.
    epsilon: Tolerance.
    iteration: Iteration number.
    gene_len: Number of unique genes among M variants.
    alpha: Prior for pcj.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.

    Returns:
    l_pdir_den: The log of the denominator of the lambda value as described in
    calculate_l_pdir.

    """

    lhs_den_norm_const = return_norm_const(gamma, pc[iteration - 1, 0, :], epsilon)
    den_density_prod = return_product_density(
        gamma, pc[iteration - 1, 0, :], pc_proposal, C
    )
    l_pdir = lhs_den_norm_const + den_density_prod
    l_pdir_gene = 0
    rhs_den_norm_const = return_norm_const(
        alpha[iteration - 1, 0], pc[iteration - 1, 0, :], epsilon
    )
    for gene_idx in range(0, gene_len):
        den_gene_density_prod = return_product_density(
            alpha[iteration - 1, 0],
            pc[iteration - 1, 0, :],
            pcj[iteration - 1, gene_idx, :],
            C,
        )
        l_pdir_gene += den_gene_density_prod + rhs_den_norm_const
    l_pdir_den = l_pdir + l_pdir_gene
    return l_pdir_den


def calculate_l_pdir(alpha, iteration, gamma, pc, pcj, epsilon, gene_len, C):

    """
    Returns the log of the transition probability of the proposal. Probability:

            /                                         J                               \
            |             /   (t - 1)           \   ____                              |
            |       p_dir \ π_0       | Γ . π_0'/ .  ||   p_dir(π_j | α . π_0')       |
            |                                      j = 1                              |
    λ = min | 1, ---------------------------------------------------------------------|
            |                                       J                                 |
            |          /               (t - 1) \   ____                     (t - 1)   |
            |    p_dir \ π_0' | Γ . π_0        / .  ||   p_dir(π_j | α . π_0        ) |
            |                                     j = 1                               | 
            \                                                                         /

    Calls helper functions for numerator and denominator.

    Parameters:
    alpha: Prior for pcj.
    iteration: Iteration number.
    gamma: Multiplicative part of the proposal value.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    epsilon: Tolerance.
    gene_len: Number of unique genes among M variants.
    C: Number of hypothesized clusters.

    Returns:
    l_pdir: Log of lambda as above.
    pc_proposal: Sampled proposal value from the Dirichlet distribution using alpha and
        previous iteration as the parameters.

    """
    
    param = alpha[iteration - 1, 0] * pc[iteration - 1, 0, :]
    if np.any(param <= 0):
        param[np.where(param <= 0)] = 1e-6
    pc_proposal = np.random.dirichlet(param)
    l_pdir_num = calculate_l_pdir_num(
        gamma, C, pc_proposal, epsilon, iteration, gene_len, alpha, pc, pcj
    )
    l_pdir_den = calculate_l_pdir_den(
        gamma, C, pc_proposal, epsilon, iteration, gene_len, alpha, pc, pcj
    )
    l_pdir = l_pdir_num - l_pdir_den
    return l_pdir, pc_proposal


def MH(
    thresh,
    accept,
    reject,
    accept_postburnin,
    reject_postburnin,
    accept_quantity,
    reject_quantity,
    iteration,
    burn,
):

    """
    Metropolis-Hastings step. With probability exp(thresh), we return the
    accept_quantity. With probability 1-exp(thresh), we return the 
    reject_quantity.

    Parameters:
    thresh: Log threshold at which we accept the new proposal.
    [accept/reject]: Tracker for acceptance rate.
    [accept/reject]_postburnin: Tracker for acceptance rate after burn-in.
    accept_quantity: Quantity to return if accepted.
    reject_quantity: Quantity to return if rejected.
    iteration: Iteration number.
    burn: Number of target burn-in iterations.

    Returns:
    [accept/reject]: Augmented tracker for acceptance rate.
    [accept/reject]_postburnin: Augmented tracker for acceptance rate after burn-in.
    quantity: One of [accept_quantity/reject_quantity] based on probability.

    """

    quantity = None
    if np.log(np.random.uniform(0, 1, size=1)[0]) < min(0, thresh):
        accept += 1
        quantity = accept_quantity
        if iteration > burn:
            accept_postburnin += 1
    else:
        reject += 1
        quantity = reject_quantity
        if iteration > burn:
            reject_postburnin += 1
    return accept, reject, accept_postburnin, reject_postburnin, quantity


def update_pc(
    alpha,
    pc,
    pcj,
    epsilon,
    gamma,
    C,
    gene_len,
    accept_mh1,
    accept_mh1_postburnin,
    reject_mh1,
    reject_mh1_postburnin,
    burn,
    iteration,
):

    """
    Calculates lambda, performs MH, and updates pc with the returned quantity.

    Parameters:
    alpha: Prior for pcj.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    epsilon: Tolerance.
    gamma: Multiplicative part of the proposal value.
    C: Number of hypothesized clusters.
    gene_len: Number of unique genes among M variants.
    [accept/reject]_mh1: Tracker for acceptance rate.
    [accept/reject]_mh1_postburnin: Tracker for acceptance rate after burn-in.
    burn: Number of target burn-in iterations.
    iteration: Iteration number.
    
    Returns:
    pc: Updated probability vector.
    [accept/reject]: Augmented tracker for acceptance rate.
    [accept/reject]_postburnin: Augmented tracker for acceptance rate after burn-in.

    """

    # Calculate lambda
    l_pdir, pc_proposal = calculate_l_pdir(
        alpha, iteration, gamma, pc, pcj, epsilon, gene_len, C
    )
    ## Metropolis-Hastings step
    accept_mh1, reject_mh1, accept_mh1_postburnin, reject_mh1_postburnin, quantity = MH(
        l_pdir,
        accept_mh1,
        reject_mh1,
        accept_mh1_postburnin,
        reject_mh1_postburnin,
        pc_proposal,
        pc[iteration - 1, 0, :],
        iteration,
        burn,
    )
    pc[iteration, 0, :] = quantity
    return pc, accept_mh1, reject_mh1, accept_mh1_postburnin, reject_mh1_postburnin


def update_pcj(alpha, pc, pcj, delta_m, gene_len, gene_vec, gene_map, iteration):

    """
    Updates per-gene probability vector dictating sharing.

                    /         / M_j           M_j              M_j         \\
                    |         | ===           ===              ===         ||
                    |         | \             \                \           ||
    π_0 ~ Dirichlet | απ_0  + | /    δ = 1,   /    δ = 2, ...  /    δ = C  ||
                    |         | ===           ===              ===         ||
                    \         \m = 1         m = 1            m = 1        //
    
    Parameters:
    alpha: Prior for pcj.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    delta_m: Indices of cluster memberships for each variant m in gene j.
    gene_len: Number of unique genes among M variants.
    gene_vec: Vector of length M of gene symbols for the M variants.
    gene_map: List of unique genes.
    iteration: Iteration number.
    
    Returns:
    pcj: Updated per-gene probability vector dictating sharing.
    
    """

    for gene_idx in range(0, gene_len):
        param_vec_shared = alpha[iteration - 1, 0] * pc[iteration, 0, :]
        for gene_iteration in range(0, len(gene_vec)):
            if gene_vec[gene_iteration] == gene_map[gene_idx]:
                param_vec_shared[int(delta_m[iteration - 1, gene_iteration])] += 1
        if np.any(param_vec_shared <= 0):
            param_vec_shared[np.where(param_vec_shared <= 0)] = 1e-6
        pcj[iteration, gene_idx, :] = np.random.dirichlet(param_vec_shared)
    return pcj


def calculate_Vjm(ses, var_idx, err_corr, Vjm_scale):

    """
    Calculates covariance structure from summary statistics that remains 
    constant across clusters.

    Parameters:
    ses: M*K matrix of standard errors.
    var_idx: Variant number in range (1 - M).
    err_corr: Correlation of errors (common vars).
    Vjm_scale: Small number added to diagonals to make sure Vjm is invertible.

    """

    dtmp = np.diag(np.array(ses[var_idx, :])[0])
    Vjm = dtmp * err_corr * dtmp + np.matlib.eye(err_corr.shape[0]) * Vjm_scale
    return Vjm


def update_delta_jm(
    betas,
    ses,
    err_corr,
    C,
    bc,
    pcj,
    delta_m,
    scales,
    maxloglkiter,
    var_idx,
    iteration,
    annot_len,
    annot_vec,
    annot_map,
    gene_len,
    gene_vec,
    gene_map,
    Vjm_scale,
):

    """
    Updates delta_jm: Indices of cluster memberships for each variant m in gene j.

               ^ 
    p_mjc = p( β_jm; bc, scales) * π_jc
    δ ~ Discrete(p_mjc)

    Parameters:
    betas: M*K matrix of effect sizes.
    ses: M*K matrix of standard errors effect sizes.
    err_corr: Correlation of errors (common vars).
    C: Number of hypothesized clusters.
    bc: Mean effect size per cluster.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    delta_m: Indices of cluster memberships for each variant m in gene j.
    scales: Scale parameters for annotations across clusters.
    maxloglkiter: Max log likelihood of the iteration.
    var_idx: Variant number in range (1 - M).
    iteration: Iteration number.
    annot_len: Number of unique annotations among M variants.
    annot_vec: Vector of length M of functional annotations.
    annot_map: List of unique annotations.
    gene_len: Number of unique genes among M variants.
    gene_vec: Vector of length M of gene symbols.
    gene_map: List of unique genes.
    Vjm_scale: Small number added to diagonals to make sure Vjm is invertible.

    Returns:
    annot_idx: Which index, in the set of unique annotations, this one corresponds to.
    delta_m: Updated indices of cluster memberships.
    
    """

    xk = np.arange(0, C)
    p_mjc, l_p_mjcu, uc = [0] * C, [0] * C, [0] * C
    var_annot = annot_vec[var_idx]
    annot_idx = [i for i in range(0, annot_len) if annot_map[i] == var_annot][0]
    gene_var = gene_vec[var_idx]
    Vjm = calculate_Vjm(ses, var_idx, err_corr, Vjm_scale)
    gene_id = [i for i in range(0, len(gene_map)) if gene_map[i] == gene_var][0]
    # Gives covariance matrix of variant effect on sets of phenotypes
    # (after fixed effect meta-analysis has been applied across all studies available)
    for c in range(0, C):
        llk2 = multivariate_normal.logpdf(
            betas[var_idx, :],
            np.sqrt(scales[iteration - 1, annot_idx]) * bc[iteration - 1, c, :],
            Vjm,
        ) + np.log(pcj[iteration, gene_id, c])
        if delta_m[iteration - 1, var_idx] == c:
            maxloglkiter[iteration - 1, 0] += llk2
        l_p_mjcu[c] += llk2
        # normalize uc - set to wc
    maxloglk = np.max(l_p_mjcu)
    for c in range(0, C):
        uc[c] = np.exp(l_p_mjcu[c] - maxloglk)
    for c in range(0, C):
        p_mjc[c] = uc[c] / np.sum(uc)
    if np.isnan(p_mjc[0]):
        wstmp = np.random.dirichlet(np.repeat(np.array([1]), C, axis=0))
        custm = stats.rv_discrete(name="custm", values=(xk, wstmp))
    else:
        custm = stats.rv_discrete(name="custm", values=(xk, p_mjc))
    delta_m[iteration, var_idx] = custm.rvs(size=1)[0]
    return annot_idx, delta_m


def update_bc(
    betas,
    ses,
    err_corr,
    C,
    M,
    delta_m,
    bc,
    scales,
    Theta_0_inv,
    iteration,
    annot_len,
    annot_vec,
    annot_map,
    Vjm_scale,
):

    """
    Updates bc using a Gibbs update from a Gaussian distribution, with mean
    and variance dependent on Theta_inv_0, scales, betas.

    Parameters:
    betas: M*K matrix of effect sizes.
    ses: M*K matrix of standard errors effect sizes.
    err_corr: Correlation of errors (common vars).
    C: Number of hypothesized clusters.
    bc: Mean effect size per cluster.
    M: Number of variants.
    delta_m: Indices of cluster memberships for each variant m in gene j.
    scales: Scale parameters for annotations across clusters.
    Theta_0_inv: Inverse of the prior estimate of genetic correlation 
        across traits.
    iteration: Iteration number.
    annot_len: Number of unique annotations among M variants.
    annot_vec: Vector of length M of functional annotations.
    annot_map: List of unique annotations.
    Vjm_scale:  Small number added to diagonals to make sure Vjm is invertible.

    Returns:
    bc: Updated mean effect size per cluster.

    """

    for c in range(1, C):
        count = 0
        mu_lhs = 0
        var_q = 0
        mu_rhs = 0 * betas[0, :]
        mu_rhs = mu_rhs.T
        for var_idx in range(0, M):
            if delta_m[iteration, var_idx] == c:
                count += 1
                varannot = annot_vec[var_idx]
                annot_idx = [
                    i for i in range(0, annot_len) if annot_map[i] == varannot
                ][0]
                Vjm = calculate_Vjm(ses, var_idx, err_corr, Vjm_scale)
                Vjminv = np.linalg.inv(Vjm)
                q1 = scales[iteration - 1, annot_idx] * Vjminv
                q2 = (
                    np.sqrt(scales[iteration - 1, annot_idx])
                    * Vjminv
                    * betas[var_idx, :].T
                )
                q3 = scales[iteration - 1, annot_idx] * Vjminv
                if count == 1:
                    mu_lhs, mu_rhs, var_q = q1, q2, q3
                else:
                    mu_lhs += q1
                    mu_rhs += q2
                    var_q += q3
        mu_lhs += Theta_0_inv
        var_q += Theta_0_inv
        mean_param = np.ravel(np.linalg.inv(mu_lhs) * mu_rhs)
        var_param = np.linalg.inv(var_q)
        bc[iteration, c, :] = np.random.multivariate_normal(mean_param, var_param)
    return bc


def update_sigma_2(
    betas,
    ses,
    xi_0,
    err_corr,
    M,
    delta_m,
    bc,
    scales,
    iteration,
    burn,
    annot_len,
    annot_vec,
    annot_map,
    accept_mh2,
    accept_mh2_postburnin,
    reject_mh2,
    reject_mh2_postburnin,
    Vjm_scale,
):

    """
    Updates sigma^2 using a MH sub-step with a random-walk proposal.
    Sequentially, for each annotation a, we sample a proposal value:
    
    σ'_a = |η|, where η ~ N(σ_a^(t-1), xi_0)

            /                               ____   ^                  ^        \
            |        p(σ_a' | sh_a, sc_a) .  ||  N(β_jm | σ_a' b_cjm, V_jm)    |
            |                          anno(v_jm)=a                            |
    λ = min | 1, --------------------------------------------------------------|
            |         (t-1)                 ____   ^         (t-1)       ^     |     
            |    p(σ_a      | sh_a, sc_a) .  ||  N(β_jm | σ_a     b_cjm, V_jm) |
            \                           anno(v_jm)=a                           /

    Parameters:
    betas: M*K matrix of effect sizes.
    ses: M*K matrix of standard errors effect sizes.
    err_corr: Correlation of errors (common vars).
    xi_0: Hyperparameter controlling the spread of the proposals.
    M: Number of variants.
    delta_m: Indices of cluster memberships for each variant m in gene j.
    bc: Mean effect size per cluster.
    scales: Scale parameters for annotations across clusters.
    iteration: Iteration number.
    burn: Number of target burn-in iterations.
    annot_len: Number of unique annotations among M variants.
    annot_vec: Vector of length M of functional annotations.
    annot_map: List of unique annotations.
    [accept/reject]_mh2: Tracker for acceptance rate.
    [accept/reject]_mh2_postburnin: Tracker for acceptance rate after burn-in.
    Vjm_scale: Small number added to diagonals to make sure Vjm is invertible.

    Returns:
    [accept/reject]_mh2: Augmented tracker for acceptance rate.
    [accept/reject]_mh2_postburnin: tracker for acceptance rate after burn-in.
    
    """

    # e) Update scale sigma^2 annot.
    for annot_idx in range(0, annot_len):
        scaleprop = abs(
            np.random.normal(np.sqrt(scales[iteration - 1, annot_idx]), xi_0, size=1)[0]
        )
        annotdata = annot_map[annot_idx]
        probnum1 = stats.invgamma.logpdf(np.power(scaleprop, 2), 1, scale=1)
        probdenom1 = stats.invgamma.logpdf(scales[iteration - 1, annot_idx], 1, scale=1)
        lnum2 = 0
        ldenom2 = 0
        for var_idx in range(0, M):
            if annot_vec[var_idx] == annotdata:
                Vjm = calculate_Vjm(ses, var_idx, err_corr, Vjm_scale)
                c = int(delta_m[iteration, var_idx])
                lnum2 += multivariate_normal.logpdf(
                    betas[var_idx, :], scaleprop * bc[iteration, c, :], Vjm
                )
                ldenom2 += multivariate_normal.logpdf(
                    betas[var_idx, :],
                    np.sqrt(scales[iteration - 1, annot_idx]) * bc[iteration, c, :],
                    Vjm,
                )
        ## Metropolis-Hastings step
        #if iteration % 100 == 0:
        #    print(probnum1, probdenom1, lnum2, ldenom2)
        (
            accept_mh2[annot_idx],
            reject_mh2[annot_idx],
            accept_mh2_postburnin[annot_idx],
            reject_mh2_postburnin[annot_idx],
            quantity,
        ) = MH(
            ((lnum2 + probnum1) - (probdenom1 + ldenom2)),
            accept_mh2[annot_idx],
            reject_mh2[annot_idx],
            accept_mh2_postburnin[annot_idx],
            reject_mh2_postburnin[annot_idx],
            np.power(scaleprop, 2),
            scales[iteration - 1, annot_idx],
            iteration,
            burn,
        )
        scales[iteration, annot_idx] = quantity
    return scales, accept_mh2, accept_mh2_postburnin, reject_mh2, reject_mh2_postburnin


def return_alpha_product_density(mult, proposal, previous, C):

    """
    Parameters:
    The density at point x given the normalization constant as defined in
    return_norm_const is:

                                   C               
                                 ____     (z_c - 1)
             p_dir(x|z) = D(z) .  ||   x_c         
                                 c = 1
 
    This function returns the product part of the density (after D(z)).

    Parameters:
    mult: Multiplicative factor in conditional. E.g. gamma / alpha.
    proposal: z_c.
    previous: x_c.
    C: Number of hypothesized clusters.
    
    Returns:
    gene_density_prod:

    """

    gene_density_prod = np.sum(
        [(mult * proposal - 1) * np.log(previous[i]) for i in range(0, C)]
    )
    return gene_density_prod


def calculate_l_adir_num(alpha_proposal, iteration, pc, pcj, epsilon, gene_len, C):

    """
    Calculates the log of the numerator of the lambda value as described in
    calculate_l_adir.

    Parameters:
    alpha_proposal: Sampled proposal value from the normal distribution using alpha and
        xi_alpha_0 as the mean and spread parameters.
    iteration: Iteration number.
    gamma: Multiplicative part of the proposal value.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    epsilon: Tolerance.
    gene_len: Number of unique genes among M variants.
    C: Number of hypothesized clusters.
    
    Returns:
    l_adir_num: The log of the numerator of the lambda value as described in
    calculate_l_adir.

    """

    ## Work on numerator (l_adir_num)
    # Set alpha_num, LHS of numerator
    alpha_num = -2 * np.log(alpha_proposal) - 1 / alpha_proposal
    rhs_num_norm_const = return_norm_const(alpha_proposal, pc[iteration, 0, :], epsilon)
    # Set l_adir_prop_gene, RHS of numerator
    l_adir_prop_gene = 0
    for gene_idx in range(0, gene_len):
        num_gene_density_prod = return_alpha_product_density(
            alpha_proposal, pc[iteration, 0, :], pcj[iteration - 1, gene_idx, :], C
        )
        l_adir_prop_gene += num_gene_density_prod + rhs_num_norm_const
    l_adir_num = alpha_num + l_adir_prop_gene
    return l_adir_num


def calculate_l_adir_den(alpha, iteration, pc, pcj, epsilon, gene_len, C):

    """
    Calculates the log of the denominator of the lambda value as described in
    calculate_l_adir.

    Parameters:
    alpha: Inverse-gamma prior for pcj.
    iteration: Iteration number.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    epsilon: Tolerance.
    gene_len: Number of unique genes among M variants.
    C: Number of hypothesized clusters.
    
    Returns:
    l_adir_den: The log of the denominator of the lambda value as described in
    calculate_l_adir.

    """

    ## Denominator (iteration - 1 in conditional)
    alpha_den = -2 * np.log(alpha[iteration - 1, 0]) - 1 / alpha[iteration - 1, 0]
    rhs_den_norm_const = return_norm_const(
        alpha[iteration - 1, 0], pc[iteration, 0, :], epsilon
    )
    # Set l_adir_prop_gene, RHS of denominator
    l_adir_gene = 0
    for gene_idx in range(0, gene_len):
        den_gene_density_prod = return_alpha_product_density(
            alpha[iteration - 1, 0],
            pc[iteration, 0, :],
            pcj[iteration - 1, gene_idx, :],
            C,
        )
        l_adir_gene += den_gene_density_prod + rhs_den_norm_const
    l_adir_den = alpha_den + l_adir_gene
    return l_adir_den


def calculate_l_adir(alpha, xi_alpha_0, iteration, pc, pcj, epsilon, gene_len, C):

    """
    Returns the log of the transition probability of the proposal. Probability:

            /                  J                                          \
            |                ____        /    (t)            (t) \        |
            |       p(α')  .  ||   p_dir \ π_j     | α' . π_0    /        |
            |               j = 1                                         |
    λ = min | 1, ---------------------------------------------------------|
            |                       J                                     |
            |      /  (t - 1) \   ____       /    (t)     (t-1)      (t)\ |
            |    p \ α        / .  ||  p_dir \ π_j     | α      . π_0   / |
            |                    j = 1                                    | 
            \                                                             /

    Calls helper functions for numerator and denominator.

    Parameters:
    alpha: Prior for pcj.
    xi_alpha_0: Fixed value controlling the variance of the proposal distribution.
    iteration: Iteration number.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    epsilon: Tolerance.
    gene_len: Number of unique genes among M variants.
    C: Number of hypothesized clusters.

    Returns:
    l_adir: Log of lambda as above.
    alpha_proposal: Sampled proposal value from the normal distribution using alpha and
        xi_alpha_0 as the mean and spread parameters.

    """

    ### Calculate acceptance probability (l_adir)
    alpha_proposal = abs(
        np.random.normal(alpha[iteration - 1, 0], xi_alpha_0, size=1)[0]
    )
    l_adir_num = calculate_l_adir_num(
        alpha_proposal, iteration, pc, pcj, epsilon, gene_len, C
    )
    l_adir_den = calculate_l_adir_den(alpha, iteration, pc, pcj, epsilon, gene_len, C)
    l_adir = l_adir_num - l_adir_den
    return l_adir, alpha_proposal


def update_alpha(
    alpha,
    pc,
    pcj,
    epsilon,
    C,
    gene_len,
    xi_alpha_0,
    iteration,
    burn,
    accept_mh3,
    accept_mh3_postburnin,
    reject_mh3,
    reject_mh3_postburnin,
):

    """
    Updates alpha using a MH sub-step with a random-walk proposal.

    Parameters:
    alpha: Prior for pcj.
    pc: C-dimensional probability vector dictating sharing of clusters across genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    epsilon: Tolerance.
    C: Number of hypothesized clusters.
    gene_len: Number of unique genes among M variants.
    xi_alpha_0: Fixed value controlling the variance of the proposal distribution.
    iteration: Iteration number.
    burn: Number of target burn-in iterations.
    [accept/reject]_mh3: Tracker for acceptance rate.
    [accept/reject]_mh3_postburnin: Tracker for acceptance rate after burn-in.
    
    Returns:
    alpha: Updated prior for pcj.
    [accept/reject]_mh3: Augmented racker for acceptance rate.
    [accept/reject]_mh3_postburnin: Augmented tracker for acceptance rate after burn-in.

    """

    # f) alpha
    l_adir, alpha_proposal = calculate_l_adir(
        alpha, xi_alpha_0, iteration, pc, pcj, epsilon, gene_len, C
    )
    ## Metropolis-Hastings step
    accept_mh3, reject_mh3, accept_mh3_postburnin, reject_mh3_postburnin, quantity = MH(
        l_adir,
        accept_mh3,
        reject_mh3,
        accept_mh3_postburnin,
        reject_mh3_postburnin,
        alpha_proposal,
        alpha[iteration - 1, :],
        iteration,
        burn,
    )
    alpha[iteration, :] = quantity
    return alpha, accept_mh3, reject_mh3, accept_mh3_postburnin, reject_mh3_postburnin


def calculate_metrics(
    outpath, fout, alpha, burn, niter, thinning, maxloglkiter, gene_len, K, M, C
):

    """
    Writes alphas out to file and returns BIC/AIC for the given C value.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    alpha: Prior for pcj.
    burn: Number of target burn-in iterations.
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    thinning: MCMC thinning parameter.
    maxloglkiter: Max log likelihood of the iteration.
    gene_len: Number of unique genes among M variants.
    K: Number of phenotypes.
    M: Number of variants.
    C: Number of hypothesized clusters.
    
    Returns:
    BIC, AIC: Goodness-of-fit measures for the C input.

    """

    alphaout = open(outpath + str(fout) + "_" + str(C) + ".mcmc.alpha", "w+")
    mean = np.mean(alpha[burn + 1 : niter + 1 : thinning, 0], axis=0)
    l95ci = np.percentile(alpha[burn + 1 : niter + 1 : thinning, 0], 2.5, axis=0)
    u95ci = np.percentile(alpha[burn + 1 : niter + 1 : thinning, 0], 97.5, axis=0)
    print("Mean alpha:")
    print(mean)
    print("alpha_m50\talpha_l95\talpha_u95", file=alphaout)
    print(("%2.2f\t%2.2f\t%2.2f") % (mean, l95ci, u95ci), file=alphaout)
    alphaout.close()
    maxllkiter = np.max(maxloglkiter[burn + 1 : niter : thinning, 0])
    BIC = -2 * maxllkiter + (K + gene_len) * (C - 1) * np.log(M)
    AIC = -2 * maxllkiter + (K + gene_len) * (C - 1) * 2
    return BIC, AIC


def scaleout_write(
    outpath, fout, annot_len, scales, burn, niter, thinning, annot_map, C
):

    """
    Writes scales out to file.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    annot_len: Number of unique annotations among M variants.
    scales: Scale parameters for annotations across clusters.
    burn: Number of target burn-in iterations.
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    thinning: MCMC thinning parameter.
    annot_map: List of unique annotations.
    
    Returns:
    None.

    """

    scaleout = open(outpath + str(fout) + "_" + str(C) + ".mcmc.scale", "w+")
    print("index\tannotation\tscale_m50\tscale_l95\tscale_u95", file=scaleout)
    for annot_idx in range(0, annot_len):
        mean = np.mean(
            np.sqrt(scales[burn + 1 : niter + 1 : thinning, annot_idx]), axis=0
        )
        l95ci = np.percentile(
            np.sqrt(scales[burn + 1 : niter + 1 : thinning, annot_idx]), 2.5, axis=0
        )
        u95ci = np.percentile(
            np.sqrt(scales[burn + 1 : niter + 1 : thinning, annot_idx]), 97.5, axis=0
        )
        print(
            ("%s\t%s\t%2.2f\t%2.2f\t%2.2f")
            % (str(annot_idx), annot_map[annot_idx], mean, l95ci, u95ci),
            file=scaleout,
        )
    scaleout.close()


def tmpbc_write(outpath, fout, K, Theta_0, C):

    """
    Writes out the Theta_0 matrix to file.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    K: Number of phenotypes.
    Theta_0: Prior estimate of genetic correlation across traits; if R_phen_use is
        true, Theta_0 = R_phen. Else, it is the identity matrix.
    
    Returns:
    None.

    """

    tmpbc = open(outpath + str(fout) + "_" + str(C) + ".theta.bc", "w+")
    for i in range(0, K):
        for j in range(0, K):
            print(Theta_0[i, j], file=tmpbc, end=" ")
        print("\n", end="", file=tmpbc)
    tmpbc.close()


def mcout_write(
    outpath, fout, M, chroff_vec, annot_vec, prot_vec, gene_vec, C, delta_m, burn, niter
):

    """
    Writes out variant cluster membership posterior probabilities to file.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    M: Number of variants.
    chroff_vec: Vector of length M of CHROM:POS:REF:ALT.
    annot_vec: Vector of length M of functional annotations for the M variants.
    prot_vec: Vector of length M of HGVSp annotations.
    gene_vec: Vector of length M of gene symbols.
    C: Number of hypothesized clusters.
    delta_m: Indices of cluster memberships for each variant m in gene j.
    burn: Number of target burn-in iterations.
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    
    Returns:
    var_prob_dict: Dictionary; key = [variant, cluster]; value = probability.

    """

    mcout = open(outpath + str(fout) + "_" + str(C) + ".mcmc.posteriors", "w+")
    var_prob_dict = {}
    mcout.write("V\tmost_severe_consequence\tHGVSp\tgene_symbol\tdescription")
    for c in range(C):
        mcout.write("\tposterior_c" + str(c))
    mcout.write("\n")
    for var_idx in range(0, M):
        mcout.write(
            str(chroff_vec[var_idx])
            + "\t"
            + str(annot_vec[var_idx])
            + "\t"
            + str(prot_vec[var_idx])
            + "\t"
            + str(gene_vec[var_idx])
            + "\t"
            + str(
                str(gene_vec[var_idx]) + ":" + str(annot_vec[var_idx]) + ":" + str(prot_vec[var_idx])
            )
        )
        for c in range(0, C):
            probclustervar = np.where(delta_m[burn + 1 : niter + 1, var_idx] == c)[
                0
            ].shape[0] / (niter - burn)
            var_prob_dict[chroff_vec[var_idx], c + 1] = probclustervar
            mcout.write("\t" + str(probclustervar))
        mcout.write("\n")
    mcout.close()
    return var_prob_dict


def probout_bcout_write(outpath, fout, C, bc, delta_m, burn, niter, thinning):

    """
    Writes cluster memberships per iteration and mean + CI of effect sizes to file.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    C: Number of hypothesized clusters.
    bc: Mean effect size per cluster.
    delta_m: Indices of cluster memberships for each variant m in gene j.
    burn: Number of target burn-in iterations.
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    thinning: MCMC thinning parameter.
    
    Returns:
    None.

    """

    probout = fout + "_" + str(C) + ".mcmc.probs"
    np.savetxt(outpath + probout, delta_m, fmt="%1.2f")
    bcout = open(outpath + str(fout) + "_" + str(C) + ".mcmc.bc", "w+")
    bcout.write("cluster")
    for phenotype in phenotypes:
        print(
            ("\t%s\t%s\t%s")
            % (phenotype + "_m50", phenotype + "_l95", phenotype + "_u95"),
            end="",
            file=bcout,
        )
    bcout.write("\n")
    for c in range(0, C):
        mean = np.mean(bc[burn + 1 : niter + 1 : thinning, c, :], axis=0)
        l95ci = np.percentile(bc[burn + 1 : niter + 1 : thinning, c, :], 2.5, axis=0)
        u95ci = np.percentile(bc[burn + 1 : niter + 1 : thinning, c, :], 97.5, axis=0)
        bcout.write(str(c))
        for phenidx in range(0, mean.shape[0]):
            print(
                ("\t%2.2f\t%2.2f\t%2.2f")
                % (mean[phenidx], l95ci[phenidx], u95ci[phenidx]),
                end="",
                file=bcout,
            )
        bcout.write("\n")
    bcout.close()


def print_rejection_rates(
    accept_mh1_postburnin,
    reject_mh1_postburnin,
    accept_mh2_postburnin,
    reject_mh2_postburnin,
    accept_mh3_postburnin,
    reject_mh3_postburnin,
):

    """
    Prints rejection rates per MH step (pc, scales, alpha).
    
    Parameters:
    [accept/reject]_mh[1/2/3]: Trackers for acceptance rates at each step.
    [accept/reject]_mh[1/2/3]_postburnin: Tracker for acceptance rates after burn-in.
    
    Returns:
    None.

    """

    rejection_rate = reject_mh1_postburnin / (
        accept_mh1_postburnin + reject_mh1_postburnin
    )
    print("Rejections and acceptances, step 1:")
    print(reject_mh1_postburnin, accept_mh1_postburnin)
    logger.info(
        ("Your acceptance rate is %2.2f")
        % (rejection_rate)
    )
    print("Rejections and acceptances, step 2:")
    print(reject_mh2_postburnin, accept_mh2_postburnin)
    print("Rejections and acceptances, step 3:")
    print(reject_mh3_postburnin, accept_mh3_postburnin)


def fdr_write(outpath, fout, fdr, M, chroff_vec, var_prob_dict, C):

    """
    Writes FDR and those variants that pass FDR to file.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    fdr: Threshold for false discovery rate (default: 0.05).
    M: Number of variants.
    chroff_vec: Vector of length M of CHROM:POS:REF:ALT.
    var_prob_dict: Dictionary; key = [variant, cluster]; value = probability.
    
    Returns:
    None.

    """

    fdrout = open(outpath + str(fout) + "_" + str(C) + ".fdr", "w+")
    print(str(fdr), file=fdrout)
    var_prob_null = []
    var_fdr_id = []
    for var_idx in range(0, M):
        var_fdr_id.append(chroff_vec[var_idx])
        var_prob_null.append(var_prob_dict[chroff_vec[var_idx], 1])
    idx_sort = sorted(range(len(var_prob_null)), key=lambda k: var_prob_null[k])
    var_prob_null_sort = [var_prob_null[i] for i in idx_sort]
    var_fdr_id_sort = [var_fdr_id[i] for i in idx_sort]
    num_fdr_tmp = 0
    counter = 0
    for i in range(0, len(var_prob_null_sort)):
        counter += 1
        num_fdr_tmp += var_prob_null_sort[i]
        fdr_tmp = num_fdr_tmp / counter
        if fdr_tmp <= fdr:
            print(var_fdr_id_sort[i], file=fdrout)
    fdrout.close()


def gene_write(outpath, fout, gene_len, gene_map, pcj, burn, niter, thinning, C):

    """
    Writes out mean and 95% CI for each gene.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    gene_len: Number of unique genes among M variants.
    gene_map: List of unique genes.
    pcj: Per-gene probability vector dictating how much sharing of clusters exists
        across genes with prior alpha.
    burn: Number of target burn-in iterations.
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    thinning: MCMC thinning parameter.
    
    Returns:
    None.

    """

    genes_dict, genedatm50, genedatl95, genedatu95 = {}, {}, {}, {}
    for gene_idx in range(0, gene_len):
        genes_dict[gene_map[gene_idx]] = gene_map[gene_idx]
        genedatm50[gene_map[gene_idx]] = np.mean(
            pcj[burn + 1 : niter + 1 : thinning, gene_idx, :], axis=0
        )
        genedatl95[gene_map[gene_idx]] = np.percentile(
            pcj[burn + 1 : niter + 1 : thinning, gene_idx, :], 2.5, axis=0
        )
        genedatu95[gene_map[gene_idx]] = np.percentile(
            pcj[burn + 1 : niter + 1 : thinning, gene_idx, :], 97.5, axis=0
        )
    geneout = open(outpath + str(fout) + "_" + str(C) + ".mcmc.gene.posteriors", "w+")
    print("gene_symbol", file=geneout, end="")
    for cluster in range(C):
        print("\tc" + str(cluster) + "_m50", file=geneout, end="")
    for cluster in range(C):
        print("\tc" + str(cluster) + "_l95", file=geneout, end="")
    for cluster in range(C):
        print("\tc" + str(cluster) + "_u95", file=geneout, end="")
    geneout.write("\n")
    for genekey in genes_dict.keys():
        print(genekey, file=geneout, end="")
        for i in range(0, len(genedatm50[genekey])):
            print(("\t%2.2f") % (genedatm50[genekey][i]), file=geneout, end="")
        for i in range(0, len(genedatl95[genekey])):
            print(("\t%2.2f") % (genedatl95[genekey][i]), file=geneout, end="")
        for i in range(0, len(genedatu95[genekey])):
            print(("\t%2.2f") % (genedatu95[genekey][i]), file=geneout, end="")
        geneout.write("\n")
    geneout.close()
    return genedatm50


def prot_write(
    outpath, fout, chroff_vec, annot_vec, prot_vec, gene_vec, protind, burn, niter, M, C
):

    """
    Writes out variant cluster membership posterior probabilities to file.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    chroff_vec: Vector of length M of CHROM:POS:REF:ALT.
    annot_vec: Vector of length M of functional annotations for the M variants.
    prot_vec: Vector of length M of HGVSp annotations.
    gene_vec: Vector of length M of gene symbols.
    protind: Protective scan array.
    burn: Number of target burn-in iterations.
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    M: Number of variants.
    
    Returns:
    None.
    
    """

    protout = open(outpath + str(fout) + "_" + str(C) + ".mcmc.protective", "w+")
    for var_idx in range(0, M):
        protout.write(
            chroff_vec[var_idx]
            + "\t"
            + annot_vec[var_idx]
            + "\t"
            + prot_vec[var_idx]
            + "\t"
            + gene_vec[var_idx]
            + "\t"
            + str(
                gene_vec[var_idx] + ":" + annot_vec[var_idx] + ":" + prot_vec[var_idx]
            )
        )
        protdattmp = np.where(protind[burn + 1 : niter + 1, var_idx] == 1)[0].shape[
            0
        ] / (niter - burn)
        protout.write("\t" + str(protdattmp))
        protout.write("\n")
    protout.close()


def mrpmm(
    betas,
    ses,
    err_corr,
    annot_vec,
    gene_vec,
    prot_vec,
    chroff_vec,
    C,
    fout,
    R_phen,
    R_phen_inv,
    phenotypes,
    R_phen_use=True,
    epsilon=1e-16,
    gamma=1,
    xi_0=1,
    xi_alpha_0=1,
    fdr=0.05,
    niter=500,
    burn=100,
    thinning=1,
    verbose=True,
    outpath="/Users/mrivas/",
    targeted=False,
    maxlor=0.693,
):

    """ 
    Performs MRPMM.
  
    Parameters: 
    betas: M*K matrix of effect sizes.
    ses: M*K matrix of standard errors effect sizes.
    err_corr: Correlation of errors (common vars).
    annot_vec: Vector of length M of consequence annotations.
    gene_vec: Vector of length M of gene symbols.
    prot_vec: Vector of length M of HGVSp annotations.
    chroff_vec: Vector of length M of CHROM:POS:REF:ALT.
    C: Number of hypothesized clusters.
    fout: Prefix for output files.
    R_phen: K*K matrix of correlations of effect sizes across phenotypes.
    R_phen_inv: Inverse of R_phen.
    phenotypes: Vector of length K of phenotype IDs.
    R_phen_use: Toggles whether or not R_phen is used.
    epsilon: Default value for calculating Gamma probabilities.
    gamma: Multiplicative part of the proposal value.
    xi_0: Hyperparameter controlling the spread of the proposals.
    xi_alpha_0: Fixed value controlling the variance of the proposal distribution.
    fdr: Threshold for false discovery rate (default: 0.05).
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    burn: Number of target burn-in iterations.
    thinning: MCMC thinning parameter.
    verbose: Prints extra materials to output files.
    outpath: Folder path prefix for output files.
    targeted: Whether or not we want to perform targeted analysis.
  
    Returns: 
    [BIC, AIC, genedat]: Measures of confidence in the cluster count.
    BIC = -2*log(p(Data | theta that maximizes C, Mc)) + vc log(n), where:
        - vc is the number of parameters (K+J)*(C-1), 
        - K is the number of phenotypes, 
        - J is the number of genes, and
        - C is the number of clusters
  
    """

    (
        betas,
        ses,
        err_corr,
        K,
        M,
        gene_len,
        annot_len,
        gene_map,
        annot_map,
        alpha,
        pc,
        pcj,
        bc,
        scales,
        delta_m,
        maxloglkiter,
        Theta_0,
        Theta_0_inv,
        accept_mh1,
        accept_mh1_postburnin,
        reject_mh1,
        reject_mh1_postburnin,
        accept_mh2,
        accept_mh2_postburnin,
        reject_mh2,
        reject_mh2_postburnin,
        accept_mh3,
        accept_mh3_postburnin,
        reject_mh3,
        reject_mh3_postburnin,
    ) = initialize_MCMC(
        niter,
        R_phen,
        R_phen_inv,
        err_corr,
        betas,
        ses,
        C,
        R_phen_use,
        gene_vec,
        annot_vec,
    )

    if targeted:
        protind = np.zeros((niter + 2, M))

    # Iterations MCMC samplers
    for iteration in range(1, niter + 1):
        if iteration % 100 == 0:
            print("Iteration " + str(iteration))
        pc, accept_mh1, reject_mh1, accept_mh1_postburnin, reject_mh1_postburnin = update_pc(
            alpha,
            pc,
            pcj,
            epsilon,
            gamma,
            C,
            gene_len,
            accept_mh1,
            accept_mh1_postburnin,
            reject_mh1,
            reject_mh1_postburnin,
            burn,
            iteration,
        )
        pcj = update_pcj(
            alpha, pc, pcj, delta_m, gene_len, gene_vec, gene_map, iteration
        )

        if pcj is None:
            return [np.nan, np.nan, {}]
        else:
            for var_idx in range(0, M):
                annot_idx, delta_m = update_delta_jm(
                    betas,
                    ses,
                    err_corr,
                    C,
                    bc,
                    pcj,
                    delta_m,
                    scales,
                    maxloglkiter,
                    var_idx,
                    iteration,
                    annot_len,
                    annot_vec,
                    annot_map,
                    gene_len,
                    gene_vec,
                    gene_map,
                    0.000001,
                )
                if targeted:
                    protbool = 0
                    protadverse = 0
                    for tmptidx in range(0, K):
                        if (
                            np.sqrt(scales[iteration - 1, annot_idx])
                            * bc[iteration - 1, delta_m[iteration, var_idx], tmptidx]
                            >= maxlor
                        ):
                            protadverse = 1
                        if (
                            np.sqrt(scales[iteration - 1, annot_idx])
                            * bc[iteration - 1, delta_m[iteration, var_idx], tmptidx]
                            < -0.1
                        ):
                            protbool = 1
                    if protbool == 1 and protadverse == 0:
                        protind[iteration, var_idx] = 1
            bc = update_bc(
                betas,
                ses,
                err_corr,
                C,
                M,
                delta_m,
                bc,
                scales,
                Theta_0_inv,
                iteration,
                annot_len,
                annot_vec,
                annot_map,
                0.000001,
            )
            scales, accept_mh2, accept_mh2_postburnin, reject_mh2, reject_mh2_postburnin = update_sigma_2(
                betas,
                ses,
                xi_0,
                err_corr,
                M,
                delta_m,
                bc,
                scales,
                iteration,
                burn,
                annot_len,
                annot_vec,
                annot_map,
                accept_mh2,
                accept_mh2_postburnin,
                reject_mh2,
                reject_mh2_postburnin,
                0.000001,
            )
            alpha, accept_mh3, reject_mh3, accept_mh3_postburnin, reject_mh3_postburnin = update_alpha(
                alpha,
                pc,
                pcj,
                epsilon,
                C,
                gene_len,
                xi_alpha_0,
                iteration,
                burn,
                accept_mh3,
                accept_mh3_postburnin,
                reject_mh3,
                reject_mh3_postburnin,
            )
    var_prob_dict = mcout_write(
        outpath,
        fout,
        M,
        chroff_vec,
        annot_vec,
        prot_vec,
        gene_vec,
        C,
        delta_m,
        burn,
        niter
    )
    if targeted:
        prot_write(
            outpath,
            fout,
            chroff_vec,
            annot_vec,
            prot_vec,
            gene_vec,
            protind,
            burn,
            niter,
            M,
            C,
        )
    fdr_write(outpath, fout, fdr, M, chroff_vec, var_prob_dict, C)
    print_rejection_rates(
        accept_mh1_postburnin,
        reject_mh1_postburnin,
        accept_mh2_postburnin,
        reject_mh2_postburnin,
        accept_mh3_postburnin,
        reject_mh3_postburnin,
    )
    if verbose:
        probout_bcout_write(outpath, fout, C, bc, delta_m, burn, niter, thinning)
        scaleout_write(
            outpath, fout, annot_len, scales, burn, niter, thinning, annot_map, C
        )
        tmpbc_write(outpath, fout, K, Theta_0, C)
        print("Mean pcj:", np.mean(pcj[burn + 1 : niter + 1 : thinning, :], axis=0))
    BIC, AIC = calculate_metrics(
        outpath, fout, alpha, burn, niter, thinning, maxloglkiter, gene_len, K, M, C
    )
    genedatm50 = gene_write(
        outpath, fout, gene_len, gene_map, pcj, burn, niter, thinning, C
    )
    return [BIC, AIC, genedatm50]


def get_betas(df, pheno1, pheno2, mode):

    """
    Retrieves betas from a pair of phenotypes using non-significant, 
        non-missing variants.
  
    Parameters: 
    df: Merged dataframe containing summary statistics.
    pheno1: First phenotype.
    pheno2: Second phenotype.
    mode: One of "null", "sig". Determines whether we want to sample from null or 
        significant variants. Useful for building out correlations of errors and 
        phenotypes respectively.
  
    Returns: 
    beta1: List of effect sizes from the first phenotype; used to compute 
        correlation.
    beta2: List of effect sizes from the second phenotype; used to compute
        correlation.
  
    """

    if ("P_" + pheno1 not in df.columns) or ("P_" + pheno2 not in df.columns):
        return [], []
    if mode == "null":
        df = df[
            (df["P_" + pheno1].astype(float) >= 1e-2)
            & (df["P_" + pheno2].astype(float) >= 1e-2)
        ]
    elif mode == "sig":
        df = df[
            (
                (df["P_" + pheno1].astype(float) <= 1e-5)
                | (df["P_" + pheno2].astype(float) <= 1e-5)
            )
        ]
    beta1 = list(df["BETA_" + pheno1])
    beta2 = list(df["BETA_" + pheno2])
    return beta1, beta2


def calculate_err(a, b, pheno1, pheno2, err_corr, err_df):

    """
    Calculates a single entry in the err_corr matrix.
    
    Parameters:
    a, b: Positional parameters within the err_corr matrix.
    pheno1: Name of first phenotype.
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
        err_beta1, err_beta2 = get_betas(err_df, pheno1, pheno2, "null")
        return pearsonr(err_beta1, err_beta2)[0] if err_beta1 else 0


def filter_for_err_corr(df):

    """
    Filters the initial dataframe for the criteria used to build err_corr.

    Parameters:
    df: Merged dataframe containing all summary statistics.

    Returns:
    df: Filtered dataframe that contains null, common, LD-independent variants.

    """

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


def build_err_corr(K, phenos, df):

    """
    Builds out a matrix of correlations between all phenotypes and studies using:
        - null (i.e. synonymous or functionally uninteresting)
        - not significant (P >= 1e-2)
        - common (MAF >= 0.01)
        - LD independent
    SNPs.

    Parameters:
    K: Number of phenotypes.
    phenos: Unique set of phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.

    Returns:
    err_corr: (S*K x S*K) matrix of correlation of errors across studies and phenotypes 
        for null variants. Used to calculate v_beta.

    """
    err_df = filter_for_err_corr(df)
    if len(err_df) == 0:
        return np.diag(np.ones(K))
    err_corr = np.zeros((K, K))
    for a, pheno1 in enumerate(phenos):
        for b, pheno2 in enumerate(phenos):
            # Location in matrix
            err_corr[a, b] = calculate_err(a, b, pheno1, pheno2, err_corr, err_df)
    err_corr = np.nan_to_num(err_corr)
    return err_corr


def calculate_phen(a, b, pheno1, pheno2, df, phenos_to_use, phen_corr):

    """
    Calculates a single entry in the phen_corr matrix.
    
    Parameters:
    a, b: Positional parameters within the R_phen matrix.
    pheno1: Name of first phenotype.
    pheno2: Name of second phenotype.
    df: Dataframe containing significant, common, LD-independent variants.
    phenos_to_use: Indicate which phenotypes to use to build R_phen.

    Returns:
    phen_corr[a, b]: One entry in the phen_corr matrix.

    """

    # If in lower triangle, do not compute; symmetric matrix
    if a > b:
        return phen_corr[b, a]
    elif a == b:
        return 1
    else:
        # if this combination of phenos doesn't exist in the map file, then nan
        if (pheno1 in phenos_to_use) and (pheno2 in phenos_to_use):
            phen_beta1, phen_beta2 = get_betas(df, pheno1, pheno2, "sig")
            return (
                pearsonr(phen_beta1, phen_beta2)[0]
                if phen_beta1 is not None
                else np.nan
            )
        else:
            return np.nan


def build_phen_corr(K, phenos, df, phenos_to_use):

    """
    Builds out a matrix of correlations between all phenotypes and studies using:
        - significant (P < 1e-5)
        - common (MAF >= 0.01)
        - LD-independent
    SNPs.

    Parameters:
    K: Number of phenotypes.
    phenos: Unique set of phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.
    phenos_to_use: Indicate which phenotypes to use to build R_phen.

    Returns:
    phen_corr: (K*K) matrix of correlations between all phenotypes and studies 
        for significant variants. Used to calculate R_phen.

    """

    phen_corr = np.zeros((K, K))
    for a, pheno1 in enumerate(phenos):
        for b, pheno2 in enumerate(phenos):
            # Location in matrix
            phen_corr[a, b] = calculate_phen(
                a, b, pheno1, pheno2, df, phenos_to_use, phen_corr
            )
    return phen_corr


def filter_for_phen_corr(df, sumstat_data):

    """
    Filters the initial dataframe for the criteria used to build R_phen.

    Parameters:
    df: Merged dataframe containing all summary statistics.
    sumstat_data: Dataframe indicating which summary statistics to use to build R_phen.

    Returns:
    df: Filtered dataframe that contains significant, common, LD-independent variants.

    """

    files_to_use = sumstat_data[sumstat_data["R_phen"] == True]
    if len(files_to_use) == 0:
        return [], []
    phenos_to_use = list(files_to_use["pheno"])
    cols_to_keep = ["V", "maf", "ld_indep"]
    for col_type in "BETA_", "P_":
        cols_to_keep.extend([col_type + pheno for pheno in phenos_to_use])
    df = df[cols_to_keep]
    # Get only LD-independent, common variants
    df = df[(df.maf >= 0.01) & (df.ld_indep == True)]
    df = df.dropna(axis=1, how="all")
    df = df.dropna()
    return df, phenos_to_use


def build_R_phen(K, phenos, df, sumstat_data):

    """
    Builds R_phen using phen_corr (calculated using the method directly above this).

    Parameters:
    K: Number of phenotypes.
    phenos: Unique set of phenotypes to use for analysis.
    df: Merged dataframe containing all relevant summary statistics.
    sumstat_data: Input file containing summary statistic paths + pheno data.

    Returns:
    R_phen: Empirical estimates of genetic correlation across phenotypes.

    """

    if K == 1:
        return np.ones((K, K))
    df, phenos_to_use = filter_for_phen_corr(df, sumstat_data)
    if len(df) == 0:
        return np.diag(np.ones(K))
    R_phen = build_phen_corr(K, phenos, df, phenos_to_use)
    return R_phen


def return_err_and_R_phen(df, phenos, K, sumstat_file):

    """ 
    Builds a matrix of correlations of errors across studies and phenotypes,
        and correlations of phenotypes.
  
    Parameters: 
    df: Dataframe that containa summary statistics.
    phenos: Unique set of phenotypes to use for analysis.
    K: Number of phenotypes.
    sumstat_file: Input file containing summary statistic paths + pheno data.

    Returns:
    err_corr: (S*K x S*K) matrix of correlation of errors across studies and phenotypes
        for null variants. Used to calculate v_beta.
    R_phen: Empirical estimates of genetic correlation across phenotypes.
  
    """

    # Sample common variants, stuff in filter + synonymous
    err_corr = build_err_corr(K, phenos, df)
    # Faster calculations, better accounts for uncertainty in estimates
    err_corr[abs(err_corr) < 0.01] = 0
    R_phen = build_R_phen(K, phenos, df, sumstat_file)
    R_phen[abs(R_phen) < 0.01] = 0
    # Get rid of any values above 0.95
    while np.max(R_phen - np.eye(len(R_phen))) > 0.9:
        R_phen = 0.9 * R_phen + 0.1 * np.diag(np.diag(R_phen))
    return err_corr, R_phen


def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def initialize_parser():

    """
    Parses inputs using argparse. 

    """

    parser = argparse.ArgumentParser(
        description="MRPMM takes in several variables that affect how it runs.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--variants",
        type=str,
        required=True,
        dest="variants",
        help="""path to file containing list of variants to include,
         one line per variant. Has a header of "V".

         format of file:

         V
         1:69081:G:C
         1:70001:G:A
         """,
    )
    parser.add_argument(
        "--phenotypes",
        type=str,
        required=True,
        dest="phenotypes",
        help="""path to tab-separated file containing list of: 
         summary statistic file paths,
         phenotypes, and
         whether or not to use the file in R_phen generation.
       
         format:
         
         path        pheno        R_phen
         /path/to/file1   pheno1     TRUE
         /path/to/file2   pheno2     FALSE
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
         and HGVSp annotations.
       
         format:
         
         V       gene_symbol     most_severe_consequence HGVSp  
         1:69081:G:C     OR4F5   5_prime_UTR_variant     ""
        """,
    )
    parser.add_argument(
        "--out_folder",
        type=str,
        default=[],
        dest="out_folder",
        help="""folder to which output(s) will be written (default: current folder).
         if folder does not exist, it will be created.""",
    )
    parser.add_argument(
        "--C",
        type=check_positive,
        nargs="+",
        default=[1],
        dest="clusters",
        help="""what number of clusters to use. must be valid ints. can input multiple
         (default: 1).""",
    )
    parser.add_argument(
        "--se_thresh",
        type=float,
        default=100,
        dest="se_thresh",
        help="""SE threshold for variant inclusion""",
    )
    parser.add_argument(
        "--maf_thresh",
        type=float,
        default=0.01,
        dest="maf_thresh",
        help="""MAF threshold for variant inclusion""",
    )
    return parser


def merge_dfs(sumstat_files, metadata):

    """
    Performs an outer merge on all of the files that have been read in;
    Annotates with metadata and sigma values.

    Parameters:
    sumstat_files: List of dataframes that contain summary statistics.
    metadata: df containing gene symbol, HGVSp, etc.

    Returns:
    df: Dataframe that is ready for err_corr/R_phen generation and for running MRPMM.

    """

    conserved_columns = ["V", "#CHROM", "POS", "REF", "ALT", "A1"]
    outer_merge = partial(pd.merge, on=conserved_columns, how="outer")
    df = reduce(outer_merge, sumstat_files)
    df = df.merge(metadata)
    to_keep = [
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "start_lost",
        "stop_lost",
        "protein_altering_variant",
        "inframe_deletion",
        "inframe_insertion",
        "splice_region_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "missense_variant",
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
    df = df[~df["most_severe_consequence"].isin(to_filter)]
    df = df[df["most_severe_consequence"].isin(to_keep)]
    return df


def rename_columns(df, pheno):

    """ 
    Renames columns such that information on phenotype is available 
        in the resultant dataframe.
  
    Additionally checks if the header contains "LOG(OR)_SE" instead of "SE".
  
    Parameters: 
    df: Input dataframe (from summary statistics).
    pheno: The phenotype from which the current summary statistic dataframe comes from.
  
    Returns: 
    df: A df with adjusted column names, e.g., "OR_white_british_cancer1085".
  
    """

    if "LOG(OR)_SE" in df.columns:
        df.rename(columns={"LOG(OR)_SE": "SE"}, inplace=True)
    columns_to_rename = ["BETA", "SE", "P"]
    renamed_columns = [(x + "_" + pheno) for x in columns_to_rename]
    df.rename(columns=dict(zip(columns_to_rename, renamed_columns)), inplace=True)
    return df


def read_in_summary_stat(path, pheno, se_thresh):

    """
    Reads in one summary statistics file.
  
    Additionally: adds a variant identifier ("V"), renames columns, and filters on 
        SE (<= 0.5).

    Parameters: 
    path: Path to file.
    pheno: Phenotype of interest.
  
    Returns: 
    df: Dataframe with renamed columns, ready for merge.

    """

    print(path)
    df = pd.read_csv(
        path,
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
    df = rename_columns(df, pheno)
    df = df[df["SE" + "_" + pheno].notnull()]
    df = df[df["SE" + "_" + pheno].astype(float) <= se_thresh]
    # Filter out HLA region
    df = df[~((df["#CHROM"] == 6) & (df["POS"].between(25477797, 36448354)))]
    return df


if __name__ == "__main__":

    parser = initialize_parser()
    args = parser.parse_args()

    print("")
    print("Valid command line arguments. Importing required packages...")
    print("")

    import array
    import math
    import pandas as pd
    import logging
    import numpy as np
    import numpy.matlib as npm
    from scipy import stats
    from scipy.stats import multivariate_normal
    from scipy.stats import invgamma
    from scipy.stats.stats import pearsonr
    from sklearn import covariance
    from functools import partial, reduce

    # Set up basic logging
    logger = logging.getLogger("Log")

    variants = pd.read_csv(args.variants, sep='\t').drop_duplicates()
    metadata = pd.read_csv(args.metadata_path, sep='\t')
    sumstat_data = pd.read_csv(args.phenotypes, sep='\t').drop_duplicates()

    chroff_vec = list(set(variants["V"]))

    phenotypes = np.unique(sumstat_data["pheno"])
    sumstat_paths = list(sumstat_data["path"])
    sumstat_files = []

    for path, pheno in zip(sumstat_paths, phenotypes):
        sumstat = read_in_summary_stat(path, pheno, args.se_thresh)
        sumstat_files.append(sumstat)

    df = merge_dfs(sumstat_files, metadata)
    err_corr, R_phen = return_err_and_R_phen(
        df, phenotypes, len(phenotypes), sumstat_data
    )

    # Filter only for variants of interest
    df = df[df["V"].isin(chroff_vec)]
    df = df[(df["maf"].astype(float) <= args.maf_thresh) & (df["maf"].astype(float) >= 0)]
    df = df[df.ld_indep == True]
    chroff_vec = list(df["V"])
    annot_vec = list(df["most_severe_consequence"])
    gene_vec = list(df["gene_symbol"])
    genes = np.unique(gene_vec)
    prot_vec = list(df['HGVSp'])
    # prot_vec = ["hgvsp"] * len(chroff_vec)

    # for now, put 0 if missing
    betas = df[["BETA_" + pheno for pheno in phenotypes]].fillna(0).values
    ses = df[["SE_" + pheno for pheno in phenotypes]].fillna(0).values

    if args.out_folder:
        out_folder = args.out_folder
    else:
        out_folder = ""

    R_phen_inv = np.linalg.inv(R_phen)
    bics, aics, genedats, clusters, log10BFs = [], [], [], [], []
    fout = "_".join(genes) + "_" + "_".join(phenotypes)

    clusters = args.clusters
    clusters = list(set(clusters + [1]))
    for C in clusters:
        [BIC, AIC, genedat] = mrpmm(
            betas,
            ses,
            err_corr,
            annot_vec,
            gene_vec,
            prot_vec,
            chroff_vec,
            C,
            fout,
            R_phen,
            R_phen_inv,
            phenotypes,
            R_phen_use=True,
            fdr=0.05,
            niter=500,
            burn=100,
            thinning=1,
            verbose=True,
            outpath=out_folder,
        )
        bics.append(BIC)
        print("BIC: " + str(BIC))
        aics.append(AIC)
        genedats.append(genedat)
    for i, C in enumerate(clusters):
        log10BFs.append((bics[0] - bics[i])/(2 * np.log(10)))
    out_df = pd.DataFrame({'num_clusters': clusters, 'BIC': bics, 'AIC': aics, 'log10BF': log10BFs})
    out_df.to_csv(out_folder + str(fout) + ".mcmc.bic.aic", sep='\t', index=False)
