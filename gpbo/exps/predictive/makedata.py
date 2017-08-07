import numpy as np
import gpbo
from gpbo.core import objectives
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', dest='index', action='store', default=0,type=int)
parser.add_argument('-d', '--dimension', dest='dimension', action='store', default=2,type=int)
args = parser.parse_args()

D = args.dimension
lb = np.array([-1.]*D)
ub = np.array([1.]*D)
#lengthscales from 0.05 to 1.5
lengthscale = int(np.random.uniform(0.05,1.2)*100)/100.
noisevar = -20.#int(np.random.uniform(-12,0)*10)/10.
maxit = 100
#outputscale will be normalized to 1
for i in range(1):
    f, xm, truemin = objectives.genmat52ojf(D,lb,ub,A=1.,ls=lengthscale,fixs=10**noisevar)
    fname = 'eihyp{}_{}_d{}_l{}_n{}_m{}.csv'.format(args.index,i,D,int(100*lengthscale),int(-10*noisevar),maxit)
    fpath = 'scratch'

    C=gpbo.core.config.eihypdefault(f,D,maxit,10**noisevar,fpath,fname)
    #C=gpbo.core.config.directdefault(f,D,20000,10**noisevar,fpath,fname)
    out = gpbo.search(C)