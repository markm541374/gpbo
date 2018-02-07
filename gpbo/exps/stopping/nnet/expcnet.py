import numpy as np
import scipy as sp
import argparse
import gpbo
import nnetfns
#parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--index', dest='index', action='store', default=0,type=int)
#args=parser.parse_args()

#kindex = gpbo.core.GPdc.MAT52
#lengthscale=0.5
D=4
f = nnetfns.nnet_cancer
fpath='results'
fname = 'eihyp_nnet_c.csv'
C=gpbo.core.config.pesfsgamma(f,D,20,1e-6,fpath,fname,ninit=20)

C.stopfn = gpbo.optimize.nstopfn
C.stoppara['nmax']=80

out = gpbo.search(C)
