#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:48:17 2018

@author: rico
"""

import src.isa.isa_fn as isa_fn 
import pandas, numpy
from optparse import OptionParser

def main(inpfile_from_wrapper=None, args_from_wrapper=None):    

    # ----- --------------
    #       option parsing
    #       -------------- -----
    # 
    # seeds are either provided as a matrix (good to reproduce results)
    # or created by ISA 
    #
    if args_from_wrapper == None:
        usage = "usage: %prog [options] arg . 20180625T1800"
        parser = OptionParser(usage)
    
        parser.add_option('-i','--inputfile',dest='inpfile',type='string')
        parser.add_option('-o','--outputfile',dest='outfile',type='string',default='isa.csv')
        parser.add_option('--seedmatrix',dest='seedfile',type='string',default=None)
        parser.add_option('--dsame',dest='dsame',type='float',default=0.80)
        parser.add_option('--dconv',dest='dconv',type='float',default=0.975)
        parser.add_option('--nseed',dest='nseed',type='int',default=100)
        parser.add_option('--seedsparsity',dest='seedsparsity',type='int',default=0)
        parser.add_option('--maxiter',dest='maxiter',type='int',default=50)
        parser.add_option('--sgc',dest='sgc',type='int',default=0)
        parser.add_option('--sgr',dest='sgr',type='int',default=1)
        parser.add_option('--thc',dest='thc',type='string',default='1,2,3')
        parser.add_option('--thr',dest='thr',type='string',default='1,2,3')
        parser.add_option('--norm',dest='norm',type='string',default='double')
        parser.add_option('--nt',action='store_true',dest='nt',default=False)
        parser.add_option('--inputhasheaders',action='store_true',dest='header',default=False)
        parser.add_option('--inputhaslabels',action='store_true',dest='label',default=False)
        parser.add_option('--nopurge',action='store_true',dest='nopurge',default=False)
        parser.add_option('--quiet',action='store_true',dest='quiet',default=False)
        parser.add_option('--nosweep',action='store_true',dest='nosweep',default=False)
        parser.add_option('--onefile',action='store_true',dest='onefile',default=False)
        parser.add_option('--gopseudo',action='store_true',dest='gopseudo',default=False)
        (options, args) = parser.parse_args()
    else:
        options = args_from_wrapper
        options.inpfile = inpfile_from_wrapper
    
    # ----- -------------------
    #       reading data matrix
    #       ------------------- -----
    # 
    # for the moment, missing values are not handled
    #
    if options.header:
        header=0
    else:
        header=None
    if options.label:
        label=0
    else:
        label=None
    A = pandas.read_csv(options.inpfile,index_col=label,header=header)
    A = A.fillna(0)
    a = A.values
    
    print('/\n/\n/ --- ISA wrapper 2018-06-25/18:00\n/\n/')
    
    print('\n/ --- matrix : \n      '+'{:d}'.format(a.shape[0])+'x'+'{:d}'.format(a.shape[1]))
    
    # ----- ----------------
    #       handling seeding
    #       ---------------- -----
    # 
    # seeds are either provided as a matrix (good to reproduce results)
    # or created by ISA 
    #
    if options.seedfile:
        S = pandas.read_csv(options.seedfile,header=None,index_col=None)
        S = S.fillna(0)
        s = S.values
        print('\n/ --- seed matrix loaded : '+\
              '{:d}'.format(s.shape[0])+'x'+'{:d}'.format(s.shape[1]))
    else:
        s = None

    # ----- -------------------------
    #       handling threshold inputs
    #       ------------------------- -----
    # 
    # allow the user to pass a list of thresholds, or a range
    # presumably, there's a better way to handle this than strings
    #
    sthr = options.thr
    if ',' in sthr:
        sthr = [float(x) for x in sthr.split(',')]
        sthr = numpy.array(sthr)
        sthr.sort()
    elif ':' in sthr:
        sthr = [float(x) for x in sthr.split(':')]
        if len(sthr)==3:
            sthr = numpy.arange(sthr[0],sthr[2],sthr[1])
            sthr.sort()
        else:
            sthr = numpy.array([1,2,3])
    else:
         sthr = [float(sthr)]
    print('\n/ --- row thresholds : ')
    print('     ',end='')
    for x in sthr:
        print(' '+'{:.2f}'.format(x),end='',flush=True)
    print('\n/')
         
    sthc = options.thc
    if ',' in sthc:
        sthc = [float(x) for x in sthc.split(',')]
        sthc = numpy.array(sthc)
        sthc.sort()
    elif ':' in sthc:
        sthc = [float(x) for x in sthc.split(':')]
        if len(sthc)==3:
            sthc = numpy.arange(sthc[0],sthc[2],sthc[1])
            sthc.sort()
        else:
            sthc = numpy.array([1,2,3])
    else:
        sthc = [float(sthc)]
    print('/ --- column thresholds : ')
    print('     ',end='')
    for x in sthc:
        print(' '+'{:.2f}'.format(x),end='',flush=True)
    print('\n/')
    # ----- ---
    #       ISA
    #       --- -----
    # 
    rsSR, csSC, sROB, sTHR, sTHC = \
    isa_fn.isa(a,rsSD=s,\
              sgr=numpy.sign(options.sgr),\
              sgc=numpy.sign(options.sgc),\
              seedsparsity=options.seedsparsity,\
              nseed=options.nseed,\
              normalisation_method=options.norm,\
              dconverged=options.dconv,\
              dsame=options.dsame,\
              sthr=sthr,\
              sthc=sthc,\
              maxiter=options.maxiter,\
              nonthreshout=options.nt,\
              quiet=options.quiet,\
              sweep=not(options.nosweep),\
              doPurge=not(options.nopurge))
    
    # ----- -----------------
    #       non-binary scores
    #       ----------------- -----
    #
    
#    if ( options.nt == 1 ):
#        rcCC, crRR = \
#        isa.normalize(a,method=options.norm)
#        csSC2 = \
#        isa.isamultiply(crRR,rsSR)
#        tm = abs(csSC2)
#        csSC2 = csSC2/tm.max(axis=0)
#        rsSR = \
#        isa.isamultiply(rcCC,csSC)
#        tm = abs(rsSR)
#        rsSR = rsSR/tm.max(axis=0)
#        csSC = csSC2
    
    # ----- ----------------------
    #       writing output to file
    #       ---------------------- -----
    # 
    # works, but could be rethought
    # needs to make sense with phenomenal
    # 
    nf = '{:0'+str(int(numpy.ceil(numpy.log10(0.5+len(sROB)))))+'d}'
    idx = ['row_threshold','col_threshold','robustness']
    ff = options.outfile
    
    if (options.onefile):
        
        col = ['M'+nf.format(x) for x in range(len(sROB))]
        
        idx.extend(A.index)
        idx.extend(A.columns)
        
        B=pandas.DataFrame(numpy.vstack([sTHR,sTHC,sROB,rsSR,csSC]),index=idx,columns=col)
        B.to_csv(ff)
        
    else:
        
        col = ['isa/M'+nf.format(x) for x in range(len(sROB))]

        tm = pandas.DataFrame(numpy.vstack([sTHR,sTHC,sROB]),index=idx,columns=col)
        tm.to_csv(ff.replace('csv','info.tsv'),sep='\t')
        
        if options.label:
            index=A.index
        else:
            index=None
        tm = pandas.DataFrame(rsSR,index=index,columns=col)
        tm.to_csv(ff.replace('csv','rowscore.tsv'),sep='\t',float_format='%.8f')
        
        if options.header:
            index=A.columns
        else:
            index=None
        
        if options.gopseudo:
            tm = pandas.DataFrame(csSC,index=index,columns=col)
            tm = tm.reset_index()
            col.insert(0,'shift')
            tm.columns = col            
            tm.to_csv(ff.replace('csv','colscore.tsv'),index=False,sep='\t',float_format='%.8f')
        else:
            tm = pandas.DataFrame(csSC,index=index,columns=col)
            tm.to_csv(ff.replace('csv','colscore.tsv'),sep='\t',float_format='%.8f')
        
    
if __name__ == '__main__':
    main()
    
