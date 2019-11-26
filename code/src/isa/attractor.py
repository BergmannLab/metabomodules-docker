#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 17:48:10 2018

@author: rico
"""

import pandas
import numpy
import glob
import matplotlib.pyplot as pl
import sys

def main(cmdarg, nkeep=179): 
  
    def corr2_coeff(A,B):
        # Rowwise mean of input arrays & subtract from input arrays themeselves
        A_mA = A - A.mean(1)[:,None]
        B_mB = B - B.mean(1)[:,None]
    
        # Sum of squares across rows
        ssA = (A_mA**2).sum(1);
        ssB = (B_mB**2).sum(1);
    
        # Finally get corr coeff
        return numpy.dot(A_mA,B_mB.T)/numpy.sqrt(numpy.dot(ssA[:,None],ssB[None]))
    
    print(cmdarg)
    fn_modu = cmdarg
    fn_pseu = 'isa.'+cmdarg.replace('colscore','pseudospectrum')
    fn_info = cmdarg.replace('colscore','info')
    fn_infoout = 'isa.'+fn_info
    fn_comp = glob.glob(cmdarg.replace('colscore','attract*colscore'))
    dcomp = 0.5
    # sig = pandas.read_csv('ps.isa_colaus1.align.decov.pseudospectrum.tsv',header=0,index_col=0,sep='\t')
    info = pandas.read_csv(fn_info,header=0,index_col=0,sep='\t')
    modu = pandas.read_csv(fn_modu,header=0,index_col=0,sep='\t')
    nmodu = modu.shape[1]
    if nmodu<179: print('not enough modules!')
    print(fn_comp[0])
    comp = pandas.read_csv(fn_comp[0],header=0,index_col=0,sep='\t')
    vcomp = comp.values
    ncomp = vcomp.shape[1]
    corr = corr2_coeff(modu.values.T,vcomp.T).T
    q0 = numpy.argmax(abs(corr),axis=1)
    q1 = numpy.max(abs(corr),axis=1)
    q0[q1<dcomp] = -1
    comp_count = numpy.array([sum(q0==j) for j in range(nmodu)])
    divver = 10**(numpy.ceil(numpy.log10(nmodu)))
    sorter = -numpy.array(comp_count)-numpy.arange(nmodu)/divver
    se = numpy.argsort(sorter)
    print(comp_count[list(se[nkeep-10:nkeep])])
    for i in range(1,len(fn_comp)):
        print(fn_comp[i])
        vcomp = pandas.read_csv(fn_comp[i],header=0,index_col=0,sep='\t')
        corr = corr2_coeff(modu.values.T,vcomp.values.T).T
        ncomp = vcomp.shape[1]+ncomp
        q0 = numpy.argmax(abs(corr),axis=1)
        q1 = numpy.max(abs(corr),axis=1)
        q0[q1<dcomp] = -1
        comp_count = numpy.array([comp_count[j]+sum(q0==j) for j in range(nmodu)])
        divver = 10**(numpy.ceil(numpy.log10(nmodu)))
        sorter = -numpy.array(comp_count)-numpy.arange(nmodu)/divver
        se = numpy.argsort(sorter)
        print(comp_count[list(se[nkeep-10:nkeep])])
    
    modu = modu.iloc[:,se[:nkeep]]
    modu.to_csv(fn_pseu,sep='\t')
    info = info.iloc[:,se[:nkeep]]
    info = info.append(pandas.Series(comp_count[se[:nkeep]]/ncomp,index=info.columns),ignore_index=True)
    info.to_csv(fn_infoout,sep='\t')
    
if __name__ == '__main__':
    cmdarg = str(sys.argv[1])
    main(cmdarg)
