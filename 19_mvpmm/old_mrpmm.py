from __future__ import print_function
from __future__ import division
from random import shuffle

# Written by Manuel A. Rivas
# Updated 11.28.2016

from optparse import OptionParser
from collections import Counter
import array
import itertools
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


if __name__ == "__main__":
    betas = []
    ses = []
    vymat = []
    annotvec = []
    genevec = []
    protvec = []
    chroffvec = []
    clusters = 
    fout = 
    Rphen = 
    Rpheninv = 
    phenidarr = 
    mrpmm(betas,ses,vymat,annotvec,genevec,protvec,chroffvec,clusters,fout,Rphen,Rpheninv,phenidarr, Rphenuse=True, fdr=.05, niter=1000,burn=100,thinning=1,verbose=True, outpath='')