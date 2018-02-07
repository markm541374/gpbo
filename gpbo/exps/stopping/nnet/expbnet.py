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
f = nnetfns.nnet_boston
fpath='results'
fname = 'eihyp_nnet_b.csv'
C=gpbo.core.config.eihypdefault(f,D,20,1e-6,fpath,fname,nrandinit=20,kindex=gpbo.core.GPdc.SQUEXP)
C.stopfn = gpbo.optimize.nstopfn
C.stoppara['nmax']=200

out = gpbo.search(C)

#fn = lambda x:f(x,**{})[0]
#r=sp.optimize.minimize(fn,[0.,0.,0.,0.])
#print(r)