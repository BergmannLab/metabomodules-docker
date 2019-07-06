#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Rico Rueedi [rrueedi]
"""
import numpy
from scipy import stats
from scipy.stats.stats import pearsonr

def pp(q,str):
    if not(q):
        print(str)
    return None
        
def buildseeds(nr,ns=100,seedsparsity=0):
    if seedsparsity<0:
        seeds = numpy.random.rand(nr,ns)
    else:
        seeds = numpy.zeros([nr,ns])
        for i in range(ns):
            sl=numpy.random.randint(0,nr-1,seedsparsity)
            seeds[sl,i]=1
    return seeds

def standardize(rc_a,axis=0):
    rc_z = stats.zscore(rc_a,axis=axis)
    return rc_z

def shuffle(X):
    ii = numpy.random.sample(X.shape).argsort(axis=0)
    ij = numpy.tile(numpy.arange(X.shape[1]), (X.shape[0], 1))
    return X[ii,ij]

def purge(A,delta):
    if A.shape[1]>1:
        B = abs(numpy.corrcoef(A,rowvar=False))>delta
        C = numpy.tril(B)
        se = C.sum(axis=1)==1
    else:
        se = [True]
    return numpy.array(se)

def normalize(rc_EE,method='double'):
    if method=='single':
        rc_CC = standardize(rc_EE,axis=1)
        cr_RR = standardize(rc_EE,axis=0).transpose()
    elif method=='singledouble':
        rc_CC = standardize(rc_EE,axis=1)
        cr_RR = standardize(rc_CC,axis=0).transpose()
    else:
        rc_CC = standardize(standardize(rc_EE,axis=0),axis=1)
        cr_RR = standardize(standardize(rc_EE,axis=1),axis=0).transpose()
    return rc_CC, cr_RR

def signature(dsA,thd,sgd):
    nd, ns = dsA.shape

    sAmean = dsA.mean(axis=0)
    sAstd = dsA.std(axis=0)
    sUP = sAmean+thd*sAstd
    sDW = sAmean-thd*sAstd

    t1 = dsA>sUP
    t2 = dsA<sDW

    if sgd == -1:  dsZ=dsA*t2
    elif sgd == 1: dsZ=dsA*t1
    else: dsZ=dsA*(t1|t2)
    
    tm = abs(dsZ)
    dsZmax = tm.max(axis=0)
    dsZmax[dsZmax==0]=1
    dsZ=dsZ/dsZmax
    return dsZ

def isamultiply(rcA,csB):
    # check that this actually corresponds to what we want:
    # it's strange given how then robustness is defined
    csBn = csB/numpy.linalg.norm(csB,axis=0)
    rsZ = numpy.dot(rcA,csBn)
    return rsZ

def diagcorr(rcA,rcB):
    cC = [pearsonr(rcA[:,i],rcB[:,i])[0] for i in range(rcA.shape[1])]
    cC = numpy.abs(numpy.array(cC))
    return cC

def testconverged(rcA,rcB,dconverged):
    cC = [pearsonr(rcA[:,i],rcB[:,i])[0]>dconverged for i in range(rcA.shape[1])]
    return numpy.array(cC)

def robustness(crRR,rcCC,rsSR,csSC):
    rsSRn = rsSR/numpy.linalg.norm(rsSR,axis=0)
    csSCn = csSC/numpy.linalg.norm(csSC,axis=0)
    sROBR=(rsSRn*isamultiply(rcCC,csSCn)).sum(axis=0)
    sROBC=(csSCn*isamultiply(crRR,rsSRn)).sum(axis=0)
    sROB=(sROBR*sROBC)**(1/2)
    return sROB

def modisa(crRR,rcCC,rsSR,thr,thc,sgr,sgc,maxiter,dconverged,dsame,doPurge=True):
    nr, nc = rcCC.shape
    csSC_proj = isamultiply(crRR,rsSR)
    csSC_prev = signature(csSC_proj,thc,sgc)
    
    i = 0

    rsSRF = numpy.empty((nr,0))
    csSCF = numpy.empty((nc,0))
    sROBF = numpy.array([])
    
    while (i<maxiter) & (csSC_prev.shape[1]>0):
        
        i=i+1

        rsSR_proj = isamultiply(rcCC,csSC_prev)
        rsSR = signature(rsSR_proj,thr,sgr)
        

        mk = (numpy.abs(rsSR)>1e-12).any(axis=0) \
            & (~numpy.isnan(rsSR).any(axis=0))
        rsSR=rsSR[:,mk]
        csSC_prev=csSC_prev[:,mk]

        csSC_proj = isamultiply(crRR,rsSR)
        csSC = signature(csSC_proj,thc,sgc)
        
        mk = (numpy.abs(csSC)>1e-12).any(axis=0) \
            & (~numpy.isnan(csSC).any(axis=0))
        csSC=csSC[:,mk]
        rsSR=rsSR[:,mk]
        csSC_prev=csSC_prev[:,mk]
        
        isconverged = numpy.array(diagcorr(csSC_prev,csSC))>dconverged

        csSC_prev=csSC

        # move converged to final
        if isconverged.any():

            rsSR_sub = rsSR[:,isconverged]
            csSC_sub = csSC[:,isconverged]
            sROB_sub = robustness(crRR,rcCC,rsSR_sub,csSC_sub)
            
            sROBF = numpy.concatenate([sROBF,sROB_sub])

            esortbyrob = numpy.argsort(-sROBF)

            sROBF = sROBF[esortbyrob]
            rsSRF = numpy.hstack([rsSRF,rsSR_sub])[:,esortbyrob]
            csSCF = numpy.hstack([csSCF,csSC_sub])[:,esortbyrob]
            
            if doPurge:
                if (nr<nc):
                    ediff=purge(rsSRF,dsame)
                else:
                    ediff=purge(csSCF,dsame)

                sROBF=sROBF[ediff]
                rsSRF=rsSRF[:,ediff]
                csSCF=csSCF[:,ediff]
            
        # keep iterating on non-converged
        csSC_prev=csSC_prev[:,~isconverged]
        
        # remove null seeds
        mk = (numpy.abs(csSC_prev)>1e-12).any(axis=0)
        csSC_prev=csSC_prev[:,mk]

    return rsSRF, csSCF, sROBF

def modcore(crRR,rcCC,rsSR,csSC,sROB,rsSEEDS,sTHR,sTHC,floorROB,thr,thc,sgr,sgc,maxiter,dconverged,dsame,doPurge=True):
    nr, nc = crRR.shape
    rsSRO,csSCO,sROBO = \
    modisa(crRR,rcCC,rsSEEDS,thr,thc,sgr,sgc,maxiter,dconverged,dsame,doPurge)

    rsSRF = rsSR
    csSCF = csSC
    sROBF = sROB
    sTHRF = sTHR
    sTHCF = sTHC

    if (rsSRO.shape[1]>0)&(doPurge):
        erob = (sROBO>floorROB)
        rsSRO = rsSRO[:,erob]
        csSCO = csSCO[:,erob]
        sROBO = sROBO[erob]

    if (rsSRO.shape[1]>0):

        if (rsSR.shape[1]>0):

            rsSRF=numpy.hstack([rsSR,rsSRO])
            csSCF=numpy.hstack([csSC,csSCO])
            sROBF=numpy.concatenate([sROB,sROBO])
            sTHRF=numpy.concatenate([sTHR,thr*numpy.ones(sROBO.shape)])
            sTHCF=numpy.concatenate([sTHC,thc*numpy.ones(sROBO.shape)])
            
            if doPurge:
                if (nr<nc):
                    ediff=purge(rsSRF,dsame)
                else:
                    ediff=purge(csSCF,dsame)
    
                sROBF=sROBF[ediff]
                sTHRF=sTHRF[ediff]
                sTHCF=sTHCF[ediff]
                rsSRF=rsSRF[:,ediff]
                csSCF=csSCF[:,ediff]

        else:

            rsSRF = rsSRO
            csSCF = csSCO
            sROBF = sROBO
            sTHRF = thr*numpy.ones(sROBF.shape)
            sTHCF = thc*numpy.ones(sROBF.shape)

    return rsSRF, csSCF, sROBF, sTHRF, sTHCF

def isa(rcA,rsSD=None,\
              sgr=0,sgc=1,\
              seedsparsity=0,nseed=100,\
              normalisation_method='double',\
              dconverged=0.975,dsame=0.90,\
              sthr=[1,2,3],sthc=[1,2,3],\
              maxiter=20,\
              sweep=True,\
              nonthreshout=False,\
              quiet=False,\
              doPurge=True):

    nr, nc = rcA.shape
    if type(sthr) in [int,float]: sthr=[sthr]
    if type(sthc) in [int,float]: sthc=[sthc]
    
    nthr = len(sthr)
    nthc = len(sthc)

    # ----- -----------
    #       build seeds
    #       ----------- -----

    if not(type(rsSD)==numpy.ndarray):
        rsSD = \
        buildseeds(nr,ns=nseed,seedsparsity=seedsparsity)
        rsSD_shuffled = \
        buildseeds(nr,ns=min([100,nseed]),seedsparsity=seedsparsity)
    else :
        rsSD_shuffled = rsSD
            

    # ----- ----------------
    #       robustness floor
    #       ---------------- -----
    #
    # the robustness floor for each threshold pair is the max of the
    # robustness of the modules obtained by running ISA
    # on the shuffled entry matrix

    rcA_shuffled = \
    shuffle(rcA)

    rcCC_shuffled, crRR_shuffled = \
    normalize(rcA_shuffled,method=normalisation_method)

    sfloorROB = numpy.zeros([nthr,nthc])
    if not(doPurge):
        pp(quiet,'\n/ --- purging is off\n')
    pp(quiet,'\n/ --- computing robustness floor\n')
    for jthr in range(nthr):
        thr=sthr[jthr]
        for jthc in range(nthc):
            thc=sthc[jthc]
            rsSR_shf, csSC_shf, sROB = \
            modisa(\
                            crRR_shuffled,rcCC_shuffled,\
                            rsSD_shuffled,thr,thc,sgr,sgc,\
                            maxiter,dconverged,dsame,doPurge)
            if len(sROB)>0:
                sfloorROB[jthr,jthc] = numpy.max(sROB)
                print('thc '+'{:.2f}'.format(thc)+', thr '+\
                      '{:.2f}'.format(thr)+', rob floor : '+\
                      '{:d}'.format(int(numpy.floor(sfloorROB[jthr,jthc])))+' / '\
                      '{:.1f}'.format(numpy.mean(sROB)/(thc*thr)**.5))
            else:
                sfloorROB[jthr,jthc] = 0
            # if not(quiet) : print('+',end='',flush=True)
        # if not(quiet) : print('/\n',end='',flush=True)
        
    # ----- -----------
    #       initial run
    #       ----------- -----
    #
    # ISA run

    rcCC, crRR = \
    normalize(rcA,method=normalisation_method)

    rsSR = numpy.empty((nr,0))
    csSC = numpy.empty((nc,0))
    sROB = numpy.array([])
    sTHR = numpy.array([])
    sTHC = numpy.array([])
    pp(quiet,'\n/ --- getting modules\n')
    for jthr in range(nthr):
        thr=sthr[jthr]
        for jthc in range(nthc):
            thc=sthc[jthc]
            floorROB=sfloorROB[jthr,jthc]
            rsSR, csSC, sROB, sTHR, sTHC = \
            modcore(\
                             crRR,rcCC,rsSR,csSC,\
                             sROB,rsSD,sTHR,sTHC,floorROB,\
                             thr,thc,sgr,sgc,maxiter,dconverged,dsame,doPurge)
            if not(quiet) : print('+',end='',flush=True)
        if not(quiet) : print('/\n',end='',flush=True)

    # ----- -----
    #       sweep
    #       ----- -----
    #
    # Use the obtained modules as seeds for a second run over all thresholds
    # Can be turned off with 'sweep' setting

    if (sweep)&((nthr>1)|(nthc>1)):
        pp(quiet,'\n/ --- sweeping\n')
        for jthr in range(nthr):
            thr=sthr[jthr]
            for jthc in range(nthc):
                thc=sthc[jthc]
                floorROB=sfloorROB[jthr,jthc]
                rsSR, csSC, sROB, sTHR, sTHC = \
                modcore(\
                                 crRR,rcCC,rsSR,csSC,\
                                 sROB,rsSR,sTHR,sTHC,floorROB,\
                                 thr,thc,sgr,sgc,maxiter,dconverged,\
                                 dsame,doPurge)
                if not(quiet) : print('+',end='',flush=True)
            if not(quiet) : print('/\n',end='',flush=True)

    if (nonthreshout):
        csSC2 = \
        isamultiply(crRR,rsSR)
        tm = abs(csSC2)
        csSC2 = csSC2/tm.max(axis=0)
        rsSR = \
        isamultiply(rcCC,csSC)
        tm = abs(rsSR)
        rsSR = rsSR/tm.max(axis=0)
        csSC = csSC2
        
    return rsSR, csSC, sROB, sTHR, sTHC
