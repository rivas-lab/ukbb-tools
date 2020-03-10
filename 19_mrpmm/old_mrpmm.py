from __future__ import print_function
from __future__ import division
from random import shuffle

# Written by Manuel A. Rivas
# Updated 11.28.2016
import argparse
from optparse import OptionParser
from collections import Counter
import array
import itertools
import pandas as pd
import math
import sys,re
import os
import logging
from scipy.stats import binom as binomial
import numpy as np
import numpy.matlib
import time
from scipy.stats import invgamma
import sklearn
import sklearn.covariance
# Set up basic logging
logger = logging.getLogger('Log')
from scipy import stats
from scipy.stats import multivariate_normal
import random

def is_pos_def(x):
    i = 0
    x = np.matrix(x)
    if np.all(np.linalg.eigvals(x) > 0):
        return True
    else:
       return False

# return BIC -2*log(p(Data | theta that maximizes C, Mc)) + vc log(n) : vc is the number of parameters (K+J)*(C-1), K is the number of phenotypes, J is the number of genes, C is the number of clusters
def mrpmm(betas,ses,vymat,annotvec,genevec,protvec,chroffvec,clusters,fout,Rphen,Rpheninv,phenidarr, Rphenuse=True, fdr=.05, niter=1000,burn=100,thinning=1,verbose=True, outpath='/Users/mrivas/'):
    print("Running MCMC algorithm...")
    epsilon = .0000000000000001
    storephensvar = []
    S = vymat
    xi0 = 1 # hyperparameter to control spread of proposals for annotation
    xialpha0 = 1
    betas = numpy.matrix(betas)
    ses = numpy.matrix(ses)
    S = numpy.matrix(S)
    Sinv = numpy.linalg.inv(S)
    # Let k be the number of clusters, where cluster 1 is the null model cluster
    C = clusters
    maxloglkiter = np.zeros((niter+2,1))
    # Let k be the number of phenotypes
    k = betas.shape[1]
    # Let m be the number of variants
    m = betas.shape[0]

    # Initialize
    #Sigma0 for alternative clusters
    if Rphenuse:
        if is_pos_def(Rphen):
            Theta0 = Rphen
            Theta0inv = Rpheninv
        else:
            Theta0  = sklearn.covariance.shrunk_covariance(Rphen)
            Theta0inv = numpy.linalg.inv(Theta0)
    else:
        Theta0  = numpy.eye(Rphen.shape[0])
        Theta0inv = numpy.linalg.inv(Theta0)
    

    #scale matrix
    geneset = set(genevec)
    genemap = list(geneset)
    annotset = set(annotvec)
    annotlen = len(annotset)
    annotmap = list(annotset)
    scales = numpy.zeros((niter+2,annotlen))
    # store the mean trait value across the clusters for individuals that are members
    bc = numpy.zeros((niter+2,C,k))
    # store the probabilities (proportions) of cluster memberships
    pc = numpy.zeros((niter+2,1,C))
    # store the probabilities (proportions) of cluster memberships for each gene
    genenum = len(set(genevec))
    pcj = numpy.zeros((niter+2,genenum,C))
    # for each iteration keep record of the variant membership
    deltam = numpy.zeros((niter+2,m))

    ###### Why are these stored separately?
    # non-normalized probabilities for each individual variant
    uc = numpy.zeros((niter+2,m,C))
    # normalized probabilities for each individual variant
    ws = numpy.zeros((niter+2,m,C))


    # for each iteration keep record of the variant membership
    tm = numpy.zeros((niter+2,m))
    #sharing parameter
    alpha = numpy.zeros((niter+2,1))
    ks = numpy.arange(1,C+1)
    sigmainvdict = {}
    sigmadict = {}
    thetadict = {}
    thetainvdict = {}
    for clusteriter in range(2,C+1):
        sigmadict[0,clusteriter] = S
        sigmainvdict[0,clusteriter] = Sinv
        thetadict[0,clusteriter] = Theta0
        thetainvdict[0,clusteriter] = Theta0inv
    # For Metropolois Hastings sub-step : keep track of acceptance rate
    acceptmh1 = 0
    rejectmh1 = 0
    acceptmh1_postburnin = 0
    rejectmh1_postburnin = 0
    acceptmh3 = 0
    rejectmh3 = 0
    acceptmh3_postburnin = 0
    rejectmh3_postburnin = 0
    acceptmh2 = [0]*annotlen
    rejectmh2 = [0]*annotlen
    acceptmh2_postburnin = [0]*annotlen
    rejectmh2_postburnin = [0]*annotlen
    # initialize \alpha : sharing of clusters across genes
    alpha[0,:] = invgamma.rvs(1,0,1,size = 1)
    # initialize pc (proportions across all variants)
    pc[0,0,:] = np.random.dirichlet([1]*C)
    # initialize pcj (proportions for each gene j)
    for geneidx in range(0,genenum):
        pcj[0,geneidx,:] = np.random.dirichlet(alpha[0,0]*pc[0,0,:])
    bc[0,0,:] = np.array([0]*k)
    for clusteridx in range(1,C):
        bc[0,clusteridx,:] = np.random.multivariate_normal(np.array([0]*k).T,Theta0)
    for scaleidx in range(0,annotlen):
        scales[0,scaleidx] = np.power(0.2,2)
    # initialize variant membership across clusters
    deltam[0,:] = np.random.randint(0,C,m)
    # Iterations MCMC samplers
    for iter in range(1,niter+1):
        gamma = 1
        if iter % 100 == 0:
            print(iter)
        ## a) Update \pi_0 : Proposal centered around the current value, Set gamma to 1 , how to set gamma?
        ## mhstep1
        pcproposal = np.random.dirichlet(alpha[iter-1,0]*pc[iter-1,0,:])
      #  lnormDprop = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pcproposal])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pcproposal])
        lnormDprop = math.lgamma(np.sum([gamma*i for i in pcproposal])) - np.sum([math.lgamma(max(gamma*i,epsilon)) for i in pcproposal])
        # second part of density
    #    densitypropb = np.sum([(alpha[iter-1,0]*pcproposal[i] - 1)*np.log(pc[iter-1,0,i]) for i in range(0,C)])
        densitypropb = np.sum([(gamma*pcproposal[i] - 1)*np.log(pc[iter-1,0,i]) for i in range(0,C)])
        lpdirprop = lnormDprop + densitypropb
        #go through each gene
        lpdirpropgene = 0
        lnormDprop = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pcproposal])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pcproposal])
        for geneidx in range(0,genenum):
            # second part of density
            densitypropb = np.sum([(alpha[iter-1,0]*pcproposal[i] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirpropgene += densitypropb + lnormDprop
        lpdirnum = lpdirprop + lpdirpropgene
        # denominator, iteration - 1 pc
     #   lnormD = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pc[iter-1,0,:]])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter-1,0,:]])
        lnormD = math.lgamma(np.sum([gamma*i for i in pc[iter-1,0,:]])) - np.sum([math.lgamma(max(gamma*i,epsilon)) for i in pc[iter-1,0,:]])
        # second part of density
        densityb = np.sum([(gamma*pc[iter-1,0,i] - 1)*np.log(pcproposal[i]) for i in range(0,C)])
        lpdir = lnormD + densityb
        #go through each gene
        lpdirgene = 0
        lnormD = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pc[iter-1,0,:]])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter-1,0,:]])
        for geneidx in range(0,genenum):
            # second part of density
            densityb = np.sum([(alpha[iter-1,0]*pc[iter-1,0,i] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirgene += densityb + lnormD
        lpdirdenom = lpdir + lpdirgene
        lpdir = lpdirnum - lpdirdenom
        ## Metropolis-Hastings step
        if numpy.log(np.random.uniform(0,1,size = 1)[0]) < min(0, lpdir):
            acceptmh1 += 1
            pc[iter,0,:] = pcproposal
            if iter > burn:
                acceptmh1_postburnin += 1
        else:
            rejectmh1 += 1
            pc[iter,0,:] = pc[iter-1,0,:]
            if iter > burn:
                rejectmh1_postburnin += 1
        # b) For each gene j = 1, ..., J update \pi_j
        for geneidx in range(0,genenum):
            paramvecshared = alpha[iter-1,0]*pc[iter,0,:]
            for geneiter in range(0,len(genevec)):
                if genevec[geneiter] == genemap[geneidx]:
                    paramvecshared[int(deltam[iter-1,geneiter])] += 1
            pcj[iter,geneidx,:] = np.random.dirichlet(paramvecshared)
        # c) Update delta_jm
        xk = numpy.arange(0,C)
        for varidx in range(0,m):
            probmjc = [0]*C
            lprobmjcu = [0]*C
            uc = [0]*C
            varannot = annotvec[varidx]
            annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
            genevar = genevec[varidx]
            geneid = [i for i in range(0,len(genemap)) if genemap[i] == genevar][0]
            atmp = np.array(ses[varidx,:])[0]
            dtmp = numpy.matlib.eye(len(atmp))
            np.fill_diagonal(dtmp,atmp)
            Vjm = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*.000001
            # Gives covariance matrix of variant  effect on sets of phenotypes (after fixed effect meta-analysis has been applied across all studies available)
            for cidx in range(0,C):
                llk2 = multivariate_normal.logpdf(betas[varidx,:],np.sqrt(scales[iter-1,annotidx])*bc[iter-1,cidx,:],Vjm) + np.log(pcj[iter,geneid,cidx])
                if deltam[iter-1,varidx] == cidx:
                    maxloglkiter[iter-1,0] += llk2
                lprobmjcu[cidx] += llk2
                        #normalize uc - set to wc
            maxloglk = numpy.max(lprobmjcu)
            for cidx in range(0,C):
                uc[cidx] = numpy.exp(lprobmjcu[cidx] - maxloglk)
            for cidx in range(0,C):
                probmjc[cidx] = uc[cidx]/numpy.sum(uc)
            if numpy.isnan(probmjc[0]):
                wstmp = numpy.random.dirichlet(numpy.repeat(numpy.array([1]),C,axis = 0))
                custm = stats.rv_discrete(name='custm',values=(xk,wstmp))
            else:
                custm = stats.rv_discrete(name='custm',values=(xk,probmjc))
            deltam[iter,varidx] = custm.rvs(size=1)[0]
        # d) Update b_c using a Gibbs update from a Gaussian distribution
        for cidx in range(1,C):
            cnt = 0
            mucurrenttmp1 = 0
            varcurrenttmp1 = 0
            mucurrenttmp2 = 0*betas[0,:]
            mucurrenttmp2 = mucurrenttmp2.T
            for varidx in range(0,m):
                if deltam[iter,varidx] == cidx:
                    cnt += 1
                    if cnt == 1:
                        varannot = annotvec[varidx]
                        annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
                        atmp = np.array(ses[varidx,:])[0]
                        dtmp = numpy.matlib.eye(len(atmp))
                        np.fill_diagonal(dtmp,atmp)
                        Vjmtmp = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*.000001
                        Vjminvtmp = np.linalg.inv(Vjmtmp)
                        mucurrenttmp1 = scales[iter-1,annotidx]*Vjminvtmp
                        mucurrenttmp2 = np.sqrt(scales[iter-1,annotidx])*Vjminvtmp*betas[varidx,:].T
                        varcurrenttmp1 = scales[iter-1,annotidx]*Vjminvtmp
                    else:
                        varannot = annotvec[varidx]
                        annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
                        atmp = np.array(ses[varidx,:])[0]
                        dtmp = numpy.matlib.eye(len(atmp))
                        np.fill_diagonal(dtmp,atmp)
                        Vjmtmp = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*.000001
                        Vjminvtmp = np.linalg.inv(Vjmtmp)
                        mucurrenttmp1 += scales[iter-1,annotidx]*Vjminvtmp
                        mucurrenttmp2 += np.sqrt(scales[iter-1,annotidx])*Vjminvtmp*betas[varidx,:].T
                        varcurrenttmp1 += scales[iter-1,annotidx]*Vjminvtmp
            mucurrenttmp1 += Theta0inv
            varcurrenttmp1 += Theta0inv
            meanparam = np.ravel(np.linalg.inv(mucurrenttmp1)*mucurrenttmp2)
            varparam = np.linalg.inv(varcurrenttmp1)
            bc[iter,cidx,:] = np.random.multivariate_normal(meanparam,varparam)
        # e) Update scale sigma^2 annot.
        for annotidx in range(0,annotlen):
            scaleprop = abs(np.random.normal(np.sqrt(scales[iter-1,annotidx]),xi0,size = 1)[0])
            annotdata = annotmap[annotidx]
            probnum1 = stats.invgamma.logpdf(np.power(scaleprop,2),1,scale=1)
            probdenom1 = stats.invgamma.logpdf(scales[iter-1,annotidx],1,scale=1)
            lnum2 = 0
            ldenom2 = 0
            for varidx in range(0,m):
                if annotvec[varidx] == annotdata:
                    atmp = np.array(ses[varidx,:])[0]
                    dtmp = numpy.matlib.eye(len(atmp))
                    np.fill_diagonal(dtmp,atmp)
                    Vjm = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*.000001
                    cidx = int(deltam[iter,varidx])
                    lnum2 += multivariate_normal.logpdf(betas[varidx,:],scaleprop*bc[iter,cidx,:],Vjm)
                    ldenom2 += multivariate_normal.logpdf(betas[varidx,:],np.sqrt(scales[iter-1,annotidx])*bc[iter,cidx,:],Vjm)
            ## Metropolis-Hastings step
            if iter % 100 == 0:
                print(probnum1,probdenom1,lnum2,ldenom2)
            if np.log(np.random.uniform(0,1,size = 1)[0]) < min(0, (lnum2 + probnum1) - (probdenom1 + ldenom2)):
                acceptmh2[annotidx] += 1
                scales[iter,annotidx] = np.power(scaleprop,2)
                if iter > burn:
                    acceptmh2_postburnin[annotidx] += 1
            else:
                rejectmh2[annotidx] += 1
                scales[iter,annotidx] = scales[iter-1,annotidx]
                if iter > burn:
                    rejectmh2_postburnin[annotidx] += 1
        # f) alpha
        alphaprop = abs(np.random.normal(alpha[iter-1,0],xialpha0,size = 1)[0])
        alphanum = -2*np.log(alphaprop) - 1/alphaprop
        alphadenom =  -2*np.log(alpha[iter-1,0]) - 1/alpha[iter-1,0]
        alphanum2 = 0
        alphadenom2 = 0
        lnormDprop = 0
        lpdirpropgene = 0
        lnormDliter = 0
        lpdirlgene = 0
        densitypropa = 0
        densitya = 0
        lnormDprop = math.lgamma(np.sum([alphaprop*i for i in pc[iter,0,:]])) - np.sum([math.lgamma(max(alphaprop*i,epsilon)) for i in pc[iter,0,:]])
        lnormDliter = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pc[iter,0,:]])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter,0,:]])
        for geneidx in range(0,genenum):
            densitypropa = np.sum([(alphaprop*pc[iter,0,:] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirpropgene += densitypropa + lnormDprop
            densitya = np.sum([(alpha[iter-1,0]*pc[iter,0,:] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirlgene += densitya + lnormDliter
        ladirnum = alphanum + lpdirpropgene
        ladirdenom = alphadenom + lpdirlgene
        ladir = ladirnum - ladirdenom
        ## Metropolis-Hastings step
        if numpy.log(np.random.uniform(0,1,size = 1)[0]) < min(0, ladir):
            acceptmh3 += 1
            alpha[iter,:] = alphaprop
            if iter > burn:
                acceptmh3_postburnin += 1
        else:
            rejectmh3 += 1
            alpha[iter,:] = alpha[iter-1,:]
            if iter > burn:
                rejectmh3_postburnin += 1
    ## Write output for input files
    mcout = open(outpath + str(fout) + '.mcmc.posteriors','w+')
    varprobdict = {}
    for varidx in range(0,m):
        mcout.write(chroffvec[varidx] + '\t' + annotvec[varidx] + '\t' + protvec[varidx] + '\t' + genevec[varidx] + '\t' + str(genevec[varidx] + ':' + annotvec[varidx] + ':' +  protvec[varidx]))
        for cidx in range(0,C):
            probclustervar = numpy.where(deltam[burn+1:niter+1,varidx] == cidx)[0].shape[0]/(niter - burn)
            varprobdict[chroffvec[varidx],cidx + 1] = probclustervar
            mcout.write('\t'  + str(probclustervar))
        mcout.write('\n')
    mcout.close()
    fdrout = open(outpath + str(fout) + '.fdr','w+')
    print(str(fdr),file = fdrout)
    varprobnull = []
    varfdrid = []
    for varidx in range(0,m):
        varfdrid.append(chroffvec[varidx])
        varprobnull.append(varprobdict[chroffvec[varidx],1])
    idxsort = sorted(range(len(varprobnull)), key=lambda k: varprobnull[k])
    varprobnullsort = [varprobnull[i] for i in idxsort]
    varfdridsort = [varfdrid[i] for i in idxsort]
    numfdrtmp = 0
    counter = 0
    varlfdr = []
    for i in range(0,len(varprobnullsort)):
        counter += 1
        numfdrtmp += varprobnullsort[i]
        fdrtmp = numfdrtmp/counter
        if fdrtmp <= fdr:
            print(varfdridsort[i], file = fdrout)
    fdrout.close()
    rejectionrate = rejectmh1_postburnin/(acceptmh1_postburnin + rejectmh1_postburnin)
    print(rejectmh1_postburnin,acceptmh1_postburnin)
    logger.info(("Your acceptance rate is %2.2f")  % ( rejectmh1_postburnin/(acceptmh1_postburnin + rejectmh1_postburnin)))
    print(rejectmh2_postburnin,acceptmh2_postburnin)
    print(rejectmh3_postburnin,acceptmh3_postburnin)
    genedatm50 = {}
    genedatl95 = {}
    genedatu95 = {}
    if verbose:
        probout = fout + '.mcmc.probs'
        numpy.savetxt(outpath + probout, deltam, fmt='%1.3f')
        bcout = open(outpath + str(fout) +  '.mcmc.bc','w+')
        bcout.write('cluster')
        for i in range(0,len(phenidarr)):
            print(("\t%s\t%s\t%s") % (phenidarr[i] + 'm50',phenidarr[i] + 'l95', phenidarr[i] + 'u95'), end = '', file = bcout)
        bcout.write('\n')
        for cidx in range(0,C):
            mean = numpy.mean(bc[burn+1:niter+1:thinning,cidx,:],axis = 0)
            l95ci = numpy.percentile(bc[burn+1:niter+1:thinning,cidx,:],2.5, axis = 0)
            u95ci = numpy.percentile(bc[burn+1:niter+1:thinning,cidx,:],97.5, axis = 0)
            bcout.write(str(cidx))
            for phenidx in range(0,mean.shape[0]):
                print(("\t%2.2f\t%2.2f\t%2.2f") % (mean[phenidx], l95ci[phenidx], u95ci[phenidx]), end = '', file = bcout)
            bcout.write('\n')
        bcout.close()
        scaleout = open(outpath + str(fout) + '.mcmc.scale','w+')
        for annotidx in range(0,annotlen):
            mean = numpy.mean(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),axis = 0)
            l95ci = numpy.percentile(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),2.5, axis = 0)
            u95ci = numpy.percentile(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),97.5, axis = 0)
            print(("%s\t%s\t%2.2f\t%2.2f\t%2.2f") % (str(annotidx),annotmap[annotidx],mean,l95ci,u95ci), file = scaleout)
        scaleout.close()
        tmpbc = open(outpath + str(fout)  + '.theta.bc', 'w+')
        for jidx in range(0,k):
            for kidx in range(0,k):
                print(Theta0[jidx,kidx], file = tmpbc,end = ' ')
            print('\n',end='',file=tmpbc)
        tmpbc.close()
        pc[0,0,:]
        print('geneset', np.mean(pcj[burn+1:niter+1:thinning,:],axis=0))
        # initialize pcj (proportions for each gene j)
        genesdict = {}
        for geneidx in range(0,genenum):
            genesdict[genemap[geneidx]] = genemap[geneidx]
            genedatm50[genemap[geneidx]] = np.mean(pcj[burn+1:niter+1:thinning,geneidx,:],axis=0)
            genedatl95[genemap[geneidx]] = np.percentile(pcj[burn+1:niter+1:thinning,geneidx,:], 2.5, axis=0)
            genedatu95[genemap[geneidx]] = np.percentile(pcj[burn+1:niter+1:thinning,geneidx,:], 97.5, axis=0)
    alphaout = open(outpath + str(fout) + '.mcmc.alpha','w+')
    mean = numpy.mean(alpha[burn+1:niter+1:thinning,0],axis = 0)
    l95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],2.5, axis = 0)
    u95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],97.5, axis = 0)
    print(mean)
    print(("%2.2f\t%2.2f\t%2.2f") % (mean,l95ci,u95ci), file = alphaout)
    alphaout.close()
    maxllkiter = np.max(maxloglkiter[burn+1:niter:thinning,0])
    BIC = -2*maxllkiter + (k+ genenum)*(C-1)*np.log(m)
    AIC = -2*maxllkiter + (k+ genenum)*(C-1)*2
    geneout = open(outpath + str(fout) + '.mcmc.gene.posteriors','w+')
    for genekey in genesdict.keys():
        print(genekey, file = geneout, end = '')
        for i in range(0,len(genedatm50[genekey])):
            print(("\t%2.2f") % (genedatm50[genekey][i]), file = geneout, end = '')
        for i in range(0,len(genedatl95[genekey])):
            print(("\t%2.2f") % (genedatl95[genekey][i]), file = geneout, end = '')
        for i in range(0,len(genedatu95[genekey])):
            print(("\t%2.2f") % (genedatu95[genekey][i]), file = geneout, end = '')
        geneout.write("\n")
    geneout.close()
   # print((k+genenum)*(C-1)*np.log(m), k,genenum,C,m,np.log(m), maxllkiter)
   # print(maxloglkiter[burn+1:niter:thinning,0])
    return [BIC,AIC,genedatm50]

def targeted(betas,ses,vymat,annotvec,genevec,protvec,chroffvec,clusters,fout,Rphen,Rpheninv,Rphenuse=True,niter=1000,burn=100,thinning=1,verbose=True, maxlor = 0.693):
    print("Running MCMC algorithm...")
    epsilon = .0000000000000001
    storephensvar = []
    S = vymat
    xi0 = 1 # hyperparameter to control spread of proposals for annotation
    xialpha0 = 1
    betas = numpy.matrix(betas)
    ses = numpy.matrix(ses)
    S = numpy.matrix(S)
    Sinv = numpy.linalg.inv(S)
    # Let k be the number of clusters, where cluster 1 is the null model cluster
    C = clusters
    maxloglkiter = np.zeros((niter+2,1))
    # Let k be the number of phenotypes
    k = betas.shape[1]
    # Let m be the number of variants
    m = betas.shape[0]
    # Initialize
    #Sigma0 for alternative clusters
    if Rphenuse:
        if is_pos_def(Rphen):
            Theta0 = Rphen
            Theta0inv = Rpheninv
        else:
            Theta0  = sklearn.covariance.shrunk_covariance(Rphen)
            Theta0inv = numpy.linalg.inv(Theta0)
    else:
        Theta0  = numpy.eye(Rphen.shape[0])
        Theta0inv = numpy.linalg.inv(Theta0)
    #scale matrix
    geneset = set(genevec)
    genemap = list(geneset)
    annotset = set(annotvec)
    annotlen = len(annotset)
    annotmap = list(annotset)
    scales = numpy.zeros((niter+2,annotlen))
    # store the mean trait value across the clusters for individuals that are members
    bc = numpy.zeros((niter+2,C,k))
    # store the probabilities (proportions) of cluster memberships
    pc = numpy.zeros((niter+2,1,C))
    # store the probabilities (proportions) of cluster memberships for each gene
    genenum = len(set(genevec))
    pcj = numpy.zeros((niter+2,genenum,C))
    # for each iteration keep record of the variant membership
    deltam = numpy.zeros((niter+2,m))
    # non-normalized probabilities for each individual variant
    uc = numpy.zeros((niter+2,m,C))
    # normalized probabilities for each individual variant
    ws = numpy.zeros((niter+2,m,C))
    # for each iteration keep record of the variant membership
    tm = numpy.zeros((niter+2,m))
    #sharing parameter
    alpha = numpy.zeros((niter+2,1))
    ks = numpy.arange(1,C+1)
    # prot scan array
    protind = numpy.zeros((niter+2,m))
    sigmainvdict = {}
    sigmadict = {}
    thetadict = {}
    thetainvdict = {}
    for clusteriter in range(2,C+1):
        sigmadict[0,clusteriter] = S
        sigmainvdict[0,clusteriter] = Sinv
        thetadict[0,clusteriter] = Theta0
        thetainvdict[0,clusteriter] = Theta0inv
    # For Metropolois Hastings sub-step : keep track of acceptance rate
    acceptmh1 = 0
    rejectmh1 = 0
    acceptmh1_postburnin = 0
    rejectmh1_postburnin = 0
    acceptmh3 = 0
    rejectmh3 = 0
    acceptmh3_postburnin = 0
    rejectmh3_postburnin = 0
    acceptmh2 = [0]*annotlen
    rejectmh2 = [0]*annotlen
    acceptmh2_postburnin = [0]*annotlen
    rejectmh2_postburnin = [0]*annotlen
    # initialize \alpha : sharing of clusters across genes
    alpha[0,:] = invgamma.rvs(1,0,1,size = 1)
    # initialize pc (proportions across all variants)
    pc[0,0,:] = np.random.dirichlet([1]*C)
    # initialize pcj (proportions for each gene j)
    for geneidx in range(0,genenum):
        pcj[0,geneidx,:] = np.random.dirichlet(alpha[0,0]*pc[0,0,:])
    bc[0,0,:] = np.array([0]*k)
    for clusteridx in range(1,C):
        bc[0,clusteridx,:] = np.random.multivariate_normal(np.array([0]*k).T,Theta0)
    for scaleidx in range(0,annotlen):
        scales[0,scaleidx] = np.power(0.2,2)
    # initialize variant membership across clusters
    deltam[0,:] = np.random.randint(0,C,m)
    # Iterations MCMC samplers
    for iter in range(1,niter+1):
        gamma = 1
        if iter % 100 == 0:
            print(iter)
        ## a) Update \pi_0 : Proposal centred around the current value, Set gamma to 1 , how to set gamma?
        ## mhstep1
        pcproposal = np.random.dirichlet(alpha[iter-1,0]*pc[iter-1,0,:])
      #  lnormDprop = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pcproposal])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pcproposal])
        lnormDprop = math.lgamma(np.sum([gamma*i for i in pcproposal])) - np.sum([math.lgamma(max(gamma*i,epsilon)) for i in pcproposal])
        # second part of density
    #    densitypropb = np.sum([(alpha[iter-1,0]*pcproposal[i] - 1)*np.log(pc[iter-1,0,i]) for i in range(0,C)])
        densitypropb = np.sum([(gamma*pcproposal[i] - 1)*np.log(pc[iter-1,0,i]) for i in range(0,C)])
        lpdirprop = lnormDprop + densitypropb
        #go through each gene
        lpdirpropgene = 0
        lnormDprop = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pcproposal])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pcproposal])
        for geneidx in range(0,genenum):
            # second part of density
            densitypropb = np.sum([(alpha[iter-1,0]*pcproposal[i] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirpropgene += densitypropb + lnormDprop
        lpdirnum = lpdirprop + lpdirpropgene
        # denominator, iteration - 1 pc
     #   lnormD = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pc[iter-1,0,:]])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter-1,0,:]])
        lnormD = math.lgamma(np.sum([gamma*i for i in pc[iter-1,0,:]])) - np.sum([math.lgamma(max(gamma*i,epsilon)) for i in pc[iter-1,0,:]])
        # second part of density
        densityb = np.sum([(gamma*pc[iter-1,0,i] - 1)*np.log(pcproposal[i]) for i in range(0,C)])
        lpdir = lnormD + densityb
        #go through each gene
        lpdirgene = 0
        lnormD = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pc[iter-1,0,:]])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter-1,0,:]])
        for geneidx in range(0,genenum):
            # second part of density
            densityb = np.sum([(alpha[iter-1,0]*pc[iter-1,0,i] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirgene += densityb + lnormD
        lpdirdenom = lpdir + lpdirgene
        lpdir = lpdirnum - lpdirdenom
        ## Metropolis-Hastings step
        if numpy.log(np.random.uniform(0,1,size = 1)[0]) < min(0, lpdir):
            acceptmh1 += 1
            pc[iter,0,:] = pcproposal
            if iter > burn:
                acceptmh1_postburnin += 1
        else:
            rejectmh1 += 1
            pc[iter,0,:] = pc[iter-1,0,:]
            if iter > burn:
                rejectmh1_postburnin += 1
        # b) For each gene j = 1, ..., J update \pi_j
        for geneidx in range(0,genenum):
            paramvecshared = alpha[iter-1,0]*pc[iter,0,:]
            for geneiter in range(0,len(genevec)):
                if genevec[geneiter] == genemap[geneidx]:
                    paramvecshared[deltam[iter-1,geneiter]] += 1
            pcj[iter,geneidx,:] = np.random.dirichlet(paramvecshared)
        # c) Update delta_jm
        xk = numpy.arange(0,C)
        for varidx in range(0,m):
            probmjc = [0]*C
            lprobmjcu = [0]*C
            uc = [0]*C
            varannot = annotvec[varidx]
            annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
            genevar = genevec[varidx]
            geneid = [i for i in range(0,len(genemap)) if genemap[i] == genevar][0]
            atmp = np.array(ses[varidx,:])[0]
            dtmp = numpy.matlib.eye(len(atmp))
            np.fill_diagonal(dtmp,atmp)
            Vjm = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*np.finfo(float).eps
            # Gives covariance matrix of variant  effect on sets of phenotypes (after fixed effect meta-analysis has been applied across all studies available)
            for cidx in range(0,C):
                llk2 = multivariate_normal.logpdf(betas[varidx,:],np.sqrt(scales[iter-1,annotidx])*bc[iter-1,cidx,:],Vjm) + np.log(pcj[iter,geneid,cidx])
                if deltam[iter-1,varidx] == cidx:
                    maxloglkiter[iter-1,0] += llk2
                lprobmjcu[cidx] += llk2
                        #normalize uc - set to wc
            maxloglk = numpy.max(lprobmjcu)
            for cidx in range(0,C):
                uc[cidx] = numpy.exp(lprobmjcu[cidx] - maxloglk)
            for cidx in range(0,C):
                probmjc[cidx] = uc[cidx]/numpy.sum(uc)
            if numpy.isnan(probmjc[0]):
                wstmp = numpy.random.dirichlet(numpy.repeat(numpy.array([1]),C,axis = 0))
                custm = stats.rv_discrete(name='custm',values=(xk,wstmp))
            else:
                custm = stats.rv_discrete(name='custm',values=(xk,probmjc))
            deltam[iter,varidx] = custm.rvs(size=1)[0]
            protbool = 0
            protadverse = 0
            for tmptidx in range(0, k):
                if np.sqrt(scales[iter-1,annotidx])*bc[iter-1,deltam[iter,varidx],tmptidx] >= maxlor:
                    protadverse = 1
                if np.sqrt(scales[iter-1,annotidx])*bc[iter-1,deltam[iter,varidx],tmptidx] < -.1:
                    protbool = 1
            if protbool == 1 and protadverse == 0:
                protind[iter,varidx] = 1
        # d) Update b_c using a Gibbs update from a Gaussian distribution
        for cidx in range(1,C):
            cnt = 0
            mucurrenttmp1 = 0
            varcurrenttmp1 = 0
            mucurrenttmp2 = 0*betas[0,:]
            mucurrenttmp2 = mucurrenttmp2.T
            for varidx in range(0,m):
                if deltam[iter,varidx] == cidx:
                    cnt += 1
                    if cnt == 1:
                        varannot = annotvec[varidx]
                        annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
                        atmp = np.array(ses[varidx,:])[0]
                        dtmp = numpy.matlib.eye(len(atmp))
                        np.fill_diagonal(dtmp,atmp)
                        Vjmtmp = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*.000001
                        Vjminvtmp = np.linalg.inv(Vjmtmp)
                        mucurrenttmp1 = scales[iter-1,annotidx]*Vjminvtmp
                        mucurrenttmp2 = np.sqrt(scales[iter-1,annotidx])*Vjminvtmp*betas[varidx,:].T
                        varcurrenttmp1 = scales[iter-1,annotidx]*Vjminvtmp
                    else:
                        varannot = annotvec[varidx]
                        annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
                        atmp = np.array(ses[varidx,:])[0]
                        dtmp = numpy.matlib.eye(len(atmp))
                        np.fill_diagonal(dtmp,atmp)
                        Vjmtmp = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*.000001
                        Vjminvtmp = np.linalg.inv(Vjmtmp)
                        mucurrenttmp1 += scales[iter-1,annotidx]*Vjminvtmp
                        mucurrenttmp2 += np.sqrt(scales[iter-1,annotidx])*Vjminvtmp*betas[varidx,:].T
                        varcurrenttmp1 += scales[iter-1,annotidx]*Vjminvtmp
            mucurrenttmp1 += Theta0inv
            varcurrenttmp1 += Theta0inv
            meanparam = np.ravel(np.linalg.inv(mucurrenttmp1)*mucurrenttmp2)
            varparam = np.linalg.inv(varcurrenttmp1)
            bc[iter,cidx,:] = np.random.multivariate_normal(meanparam,varparam)
        # e) Update scale sigma^2 annot.
        for annotidx in range(0,annotlen):
            scaleprop = abs(np.random.normal(np.sqrt(scales[iter-1,annotidx]),xi0,size = 1)[0])
            annotdata = annotmap[annotidx]
            probnum1 = stats.invgamma.logpdf(np.power(scaleprop,2),1,scale=1)
            probdenom1 = stats.invgamma.logpdf(scales[iter-1,annotidx],1,scale=1)
            lnum2 = 0
            ldenom2 = 0
            for varidx in range(0,m):
                if annotvec[varidx] == annotdata:
                    atmp = np.array(ses[varidx,:])[0]
                    dtmp = numpy.matlib.eye(len(atmp))
                    np.fill_diagonal(dtmp,atmp)
                    Vjm = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*np.finfo(float).eps
                    cidx = deltam[iter,varidx]
                    lnum2 += multivariate_normal.logpdf(betas[varidx,:],scaleprop*bc[iter,cidx,:],Vjm)
                    ldenom2 += multivariate_normal.logpdf(betas[varidx,:],np.sqrt(scales[iter-1,annotidx])*bc[iter,cidx,:],Vjm)
            ## Metropolis-Hastings step
            if iter % 100 == 0:
                print(probnum1,probdenom1,lnum2,ldenom2)
            if np.log(np.random.uniform(0,1,size = 1)[0]) < min(0, (lnum2 + probnum1) - (probdenom1 + ldenom2)):
                acceptmh2[annotidx] += 1
                scales[iter,annotidx] = np.power(scaleprop,2)
                if iter > burn:
                    acceptmh2_postburnin[annotidx] += 1
            else:
                rejectmh2[annotidx] += 1
                scales[iter,annotidx] = scales[iter-1,annotidx]
                if iter > burn:
                    rejectmh2_postburnin[annotidx] += 1
        # f) alpha
        alphaprop = abs(np.random.normal(alpha[iter-1,0],xialpha0,size = 1)[0])
        alphanum = -2*np.log(alphaprop) - 1/alphaprop
        alphadenom =  -2*np.log(alpha[iter-1,0]) - 1/alpha[iter-1,0]
        alphanum2 = 0
        alphadenom2 = 0
        lnormDprop = 0
        lpdirpropgene = 0
        lnormDliter = 0
        lpdirlgene = 0
        densitypropa = 0
        densitya = 0
        lnormDprop = math.lgamma(np.sum([alphaprop*i for i in pc[iter,0,:]])) - np.sum([math.lgamma(max(alphaprop*i,epsilon)) for i in pc[iter,0,:]])
        lnormDliter = math.lgamma(np.sum([alpha[iter-1,0]*i for i in pc[iter,0,:]])) - np.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter,0,:]])
        for geneidx in range(0,genenum):
            densitypropa = np.sum([(alphaprop*pc[iter,0,:] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirpropgene += densitypropa + lnormDprop
            densitya = np.sum([(alpha[iter-1,0]*pc[iter,0,:] - 1)*np.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirlgene += densitya + lnormDliter
        ladirnum = alphanum + lpdirpropgene
        ladirdenom = alphadenom + lpdirlgene
        ladir = ladirnum - ladirdenom
        ## Metropolis-Hastings step
        if numpy.log(np.random.uniform(0,1,size = 1)[0]) < min(0, ladir):
            acceptmh3 += 1
            alpha[iter,:] = alphaprop
            if iter > burn:
                acceptmh3_postburnin += 1
        else:
            rejectmh3 += 1
            alpha[iter,:] = alpha[iter-1,:]
            if iter > burn:
                rejectmh3_postburnin += 1
    ## Write output for input files
    mcout = open(outpath + str(fout) + '.mcmc.posteriors','w+')
    for varidx in range(0,m):
        mcout.write(chroffvec[varidx] + '\t' + annotvec[varidx] + '\t' + protvec[varidx] + '\t' + genevec[varidx] + '\t' + str(genevec[varidx] + ':' + annotvec[varidx] + ':' +  protvec[varidx]))
        for cidx in range(0,C):
            probclustervar = numpy.where(deltam[burn+1:niter+1,varidx] == cidx)[0].shape[0]/(niter - burn)
            mcout.write('\t'  + str(probclustervar))
        mcout.write('\n')
    mcout.close()
    ## Write output for input files
    protout = open(outpath + str(fout) + '.mcmc.protective','w+')
    for varidx in range(0,m):
        protout.write(chroffvec[varidx] + '\t' + annotvec[varidx] + '\t' + protvec[varidx] + '\t' + genevec[varidx] + '\t' + str(genevec[varidx] + ':' + annotvec[varidx] + ':' +  protvec[varidx]))
        protdattmp = numpy.where(protind[burn+1:niter+1,varidx] == 1)[0].shape[0]/(niter - burn)
        protout.write('\t'  + str(protdattmp))
        protout.write('\n')
    protout.close()
    rejectionrate = rejectmh1_postburnin/(acceptmh1_postburnin + rejectmh1_postburnin)
    print(rejectmh1_postburnin,acceptmh1_postburnin)
    logger.info(("Your acceptance rate is %2.2f")  % ( rejectmh1_postburnin/(acceptmh1_postburnin + rejectmh1_postburnin)))
    print(rejectmh2_postburnin,acceptmh2_postburnin)
    print(rejectmh3_postburnin,acceptmh3_postburnin)
    genedat = {}
    if verbose:
        probout = fout + '.mcmc.probs'
        numpy.savetxt(outpath + probout, deltam, fmt='%1.3f')
        bcout = open(outpath + str(fout) +  '.mcmc.bc','w+')
        for cidx in range(0,C):
            mean = numpy.mean(bc[burn+1:niter+1:thinning,cidx,:],axis = 0)
            l95ci = numpy.percentile(bc[burn+1:niter+1:thinning,cidx,:],2.5, axis = 0)
            u95ci = numpy.percentile(bc[burn+1:niter+1:thinning,cidx,:],97.5, axis = 0)
            bcout.write(str(cidx))
            for phenidx in range(0,mean.shape[0]):
                print(("\t%2.2f\t%2.2f\t%2.2f") % (mean[phenidx], l95ci[phenidx], u95ci[phenidx]), end = '', file = bcout)
            bcout.write('\n')
        bcout.close()
        scaleout = open(outpath + str(fout) + '.mcmc.scale','w+')
        for annotidx in range(0,annotlen):
            mean = numpy.mean(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),axis = 0)
            l95ci = numpy.percentile(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),2.5, axis = 0)
            u95ci = numpy.percentile(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),97.5, axis = 0)
            print(("%s\t%s\t%2.2f\t%2.2f\t%2.2f") % (str(annotidx),annotmap[annotidx],mean,l95ci,u95ci), file = scaleout)
        scaleout.close()
        tmpbc = open(outpath + str(fout)  + '.theta.bc', 'w+')
        for jidx in range(0,k):
            for kidx in range(0,k):
                print(Theta0[jidx,kidx], file = tmpbc,end = ' ')
            print('\n',end='',file=tmpbc)
        tmpbc.close()
        pc[0,0,:]
        print('geneset', np.mean(pcj[burn+1:niter+1:thinning,:],axis=0))
        # initialize pcj (proportions for each gene j)
        for geneidx in range(0,genenum):
            genedat[genemap[geneidx]] = np.mean(pcj[burn+1:niter+1:thinning,geneidx,:],axis=0)
    alphaout = open(outpath + str(fout) + '.mcmc.alpha','w+')
    mean = numpy.mean(alpha[burn+1:niter+1:thinning,0],axis = 0)
    l95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],2.5, axis = 0)
    u95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],97.5, axis = 0)
    print(mean)
    print(("%2.2f\t%2.2f\t%2.2f") % (mean,l95ci,u95ci), file = alphaout)
    alphaout.close()
    maxllkiter = np.max(maxloglkiter[burn+1:niter:thinning,0])
    BIC = -2*maxllkiter + (k+ genenum)*(C-1)*np.log(m)
    AIC = -2*maxllkiter + (k+ genenum)*(C-1)*2
   # print((k+genenum)*(C-1)*np.log(m), k,genenum,C,m,np.log(m), maxllkiter)
   # print(maxloglkiter[burn+1:niter:thinning,0])
    return [BIC,AIC,genedat]

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

def initialize_parser():

    """
    Parses inputs using argparse. 

    """

    parser = argparse.ArgumentParser(
        description="MRP takes in several variables that affect how it runs.",
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
        "--clusters",
        type=int,
        nargs="+",
        required=True,
        dest="clusters",
        help="""number of clusters hypothesized - can input more than one,
         e.g. --clusters 3 2 4. MUST be integers.""",
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
        "--fout",
        type=str,
        required=True,
        dest="fout",
        help="""file prefix for output.""",
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

def read_in_summary_stat(path, pheno):

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
    df = df[df["SE" + "_" + pheno].astype(float) <= 0.2]
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

    variants = pd.read_table(args.variants).drop_duplicates()
    metadata = pd.read_table(args.metadata_path)
    sumstat_data = pd.read_table(args.phenotypes).drop_duplicates()

    chroff_vec = list(set(variants["V"]))

    phenotypes = np.unique(sumstat_data["pheno"])
    sumstat_paths = list(sumstat_data["path"])
    sumstat_files = []

    for path, pheno in zip(sumstat_paths, phenotypes):
        sumstat = read_in_summary_stat(path, pheno)
        sumstat_files.append(sumstat)

    df = merge_dfs(sumstat_files, metadata)

    err_corr, R_phen = return_err_and_R_phen(
        df, phenotypes, len(phenotypes), sumstat_data
    )

    # Filter only for variants of interest
    df = df[df["V"].isin(chroff_vec)]
    chroff_vec = list(df["V"])
    annot_vec = list(df["most_severe_consequence"])
    gene_vec = list(df["gene_symbol"])
    # prot_vec = list(metadata['HGVSp'])
    prot_vec = ["hgvsp1", "hgvsp2", "hgvsp3", "hgvsp4"]

    # for now, put 0 if missing
    betas = df[["BETA_" + pheno for pheno in phenotypes]].fillna(0).values
    ses = df[["SE_" + pheno for pheno in phenotypes]].fillna(0).values

    if args.out_folder:
        out_folder = args.out_folder
    else:
        out_folder = ""

    R_phen_inv = np.linalg.inv(R_phen)
    bics = []
    aics = []
    genedats = []


    for C in list(set(args.clusters)):
        [BIC, AIC, genedat] = mrpmm(
            betas,
            ses,
            err_corr,
            annot_vec,
            gene_vec,
            prot_vec,
            chroff_vec,
            C,
            args.fout,
            R_phen,
            R_phen_inv,
            phenotypes,
            Rphenuse=True,
            fdr=0.05,
            niter=1000,
            burn=100,
            thinning=1,
            verbose=True,
            outpath="",
        )
    out_df = pd.DataFrame({'num_clusters': args.clusters, 'BIC': bics, 'AIC': aics})
    out_df.to_csv(out_folder + str(args.fout) + "_" + str(C) + ".mcmc.bic.aic", sep='\t', index=False)
