# coding: utf-8
from __future__ import print_function
from __future__ import division
from random import shuffle
from optparse import OptionParser
from collections import Counter
import array
import itertools
import math
import sys, re
import os
import pandas as pd
import logging
from scipy.stats import binom as binomial
import numpy as np
import numpy.matlib as npm
import time
from scipy.stats import invgamma
import sklearn
import sklearn.covariance

# Written by Manuel A. Rivas
# Updated 02.11.2020, Guhan R. Venkataraman

# Set up basic logging
logger = logging.getLogger("Log")
from scipy import stats
from scipy.stats import multivariate_normal
import random


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

    print("Running MCMC algorithm...")
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
            Theta_0 = sklearn.covariance.shrunk_covariance(R_phen)
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

    pc_proposal = np.random.dirichlet(alpha[iteration - 1, 0] * pc[iteration - 1, 0, :])
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
        if iteration % 100 == 0:
            print(probnum1, probdenom1, lnum2, ldenom2)
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
    thinning: MCMC thinning paramter.
    maxloglkiter: Max log likelihood of the iteration.
    gene_len: Number of unique genes among M variants.
    K: Number of phenotypes.
    M: Number of variants.
    C: Number of hypothesized clusters.
    
    Returns:
    BIC, AIC: Goodness-of-fit measures for the C input.

    """

    alphaout = open(outpath + str(fout) + ".mcmc.alpha", "w+")
    mean = np.mean(alpha[burn + 1 : niter + 1 : thinning, 0], axis=0)
    l95ci = np.percentile(alpha[burn + 1 : niter + 1 : thinning, 0], 2.5, axis=0)
    u95ci = np.percentile(alpha[burn + 1 : niter + 1 : thinning, 0], 97.5, axis=0)
    print(mean)
    print(("%2.2f\t%2.2f\t%2.2f") % (mean, l95ci, u95ci), file=alphaout)
    alphaout.close()
    maxllkiter = np.max(maxloglkiter[burn + 1 : niter : thinning, 0])
    BIC = -2 * maxllkiter + (K + gene_len) * (C - 1) * np.log(M)
    AIC = -2 * maxllkiter + (K + gene_len) * (C - 1) * 2
    return BIC, AIC


def scaleout_write(outpath, fout, annot_len, scales, burn, niter, thinning, annot_map):

    """
    Writes scales out to file.

    Parameters:
    outpath: Folder path prefix for output files.
    fout: Prefix for output files.
    annot_len: Number of unique annotations among M variants.
    scales: Scale parameters for annotations across clusters.
    burn: Number of target burn-in iterations.
    niter: Number of iterations for Markov Chain Monte Carlo (MCMC).
    thinning: MCMC thinning paramter.
    annot_map: List of unique annotations.
    
    Returns:
    None.

    """

    scaleout = open(outpath + str(fout) + ".mcmc.scale", "w+")
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


def tmpbc_write(outpath, fout, K, Theta_0):
    
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

    tmpbc = open(outpath + str(fout) + ".theta.bc", "w+")
    for i in range(0, K):
        for j in range(0, K):
            print(Theta_0[i, j], file=tmpbc, end=" ")
        print("\n", end="", file=tmpbc)
    tmpbc.close()

def mcout_write(outpath, fout, M, chroff_vec, annot_vec, prot_vec, gene_vec, C, delta_m, burn, niter):

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

    mcout = open(outpath + str(fout) + '.mcmc.posteriors','w+')
    var_prob_dict = {}
    for var_idx in range(0, M):
        mcout.write(
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
    thinning: MCMC thinning paramter.
    
    Returns:
    None.

    """

    probout = fout + ".mcmc.probs"
    np.savetxt(outpath + probout, delta_m, fmt="%1.3f")
    bcout = open(outpath + str(fout) + ".mcmc.bc", "w+")
    bcout.write("cluster")
    for phenotype in phenotypes:
        print(
            ("\t%s\t%s\t%s")
            % (phenotype + "m50", phenotype + "l95", phenotype + "u95"),
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

    rejectionrate = reject_mh1_postburnin / (
        accept_mh1_postburnin + reject_mh1_postburnin
    )
    print(rejectionrate)
    print(reject_mh1_postburnin, accept_mh1_postburnin)
    logger.info(
        ("Your acceptance rate is %2.2f")
        % (reject_mh1_postburnin / (accept_mh1_postburnin + reject_mh1_postburnin))
    )
    print(reject_mh2_postburnin, accept_mh2_postburnin)
    print(reject_mh3_postburnin, accept_mh3_postburnin)


def fdr_write(outpath, fout, fdr, M, chroff_vec, var_prob_dict):

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

    fdrout = open(outpath + str(fout) + ".fdr", "w+")
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


def gene_write(outpath, fout, gene_len, gene_map, pcj, burn, niter, thinning):

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
    thinning: MCMC thinning paramter.
    
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
    geneout = open(outpath + str(fout) + ".mcmc.gene.posteriors", "w+")
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
    outpath, fout, chroff_vec, annot_vec, prot_vec, gene_vec, protind, burn, niter, M
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

    protout = open(outpath + str(fout) + ".mcmc.protective", "w+")
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
    niter=1000,
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
    thinning: MCMC thinning paramter.
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
            print(iteration)
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
    var_prob_dict = mcout_write(outpath, fout, M, chroff_vec, annot_vec, prot_vec, gene_vec, C, delta_m, burn, niter)
    if targeted:
        prot_write(outpath, fout, chroff_vec, annot_vec, prot_vec, gene_vec, protind, burn, niter, M)
    fdr_write(outpath, fout, fdr, M, chroff_vec, var_prob_dict)
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
            outpath, fout, annot_len, scales, burn, niter, thinning, annot_map
        )
        tmpbc_write(outpath, fout, K, Theta_0)
        print("gene_set", np.mean(pcj[burn + 1 : niter + 1 : thinning, :], axis=0))
    BIC, AIC  = calculate_metrics(
        outpath, fout, alpha, burn, niter, thinning, maxloglkiter, gene_len, K, M, C
    )
    genedatm50 = gene_write(outpath, fout, gene_len, gene_map, pcj, burn, niter, thinning)
    return [BIC, AIC, genedatm50]


if __name__ == "__main__":
    ang = pd.read_table("ANGPTL7.tsv")
    # for now, put 0 if missing
    betas = (
        ang[
            [
                "BETA_white_british_HC276",
                "BETA_white_british_INI5255",
                "BETA_white_british_INI5257",
            ]
        ]
        .fillna(0)
        .values
    )
    ses = (
        ang[
            [
                "SE_white_british_HC276",
                "SE_white_british_INI5255",
                "SE_white_british_INI5257",
            ]
        ]
        .fillna(0)
        .values
    )
    # vymat = err_corr
    err_corr = np.array(
        [
            [1, 0.06741325, 0.03541408],
            [0.06741325, 1, 0.56616657],
            [0.03541408, 0.56616657, 1],
        ]
    )
    annot_vec = [
        "missense_variant",
        "missense_variant",
        "missense_variant",
        "stop_gained",
    ]
    gene_vec = ["ANGPTL7"] * len(annot_vec)
    prot_vec = ["hgvsp1", "hgvsp2", "hgvsp3", "hgvsp4"]
    chroff_vec = [
        "1:11252369:G:A",
        "1:11253684:G:T",
        "1:11252357:A:G",
        "1:11253688:C:T",
    ]
    C = 2
    fout = "ANGPTL7_test"
    R_phen = np.array(
        [
            [1, 0.8568072, 0.61924757],
            [0.8568072, 1, 0.82642932],
            [0.61924757, 0.82642932, 1],
        ]
    )
    R_phen_inv = np.linalg.inv(R_phen)
    phenotypes = [
        "HC276",
        "INI5255",
        "INI5257",
    ]
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
        niter=1000,
        burn=100,
        thinning=1,
        verbose=True,
        outpath="",
    )