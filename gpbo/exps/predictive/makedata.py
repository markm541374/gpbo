import numpy as np
import gpbo
import os
from gpbo.core import objectives
from gpbo.core import GPdc
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', dest='index', action='store', default=0,type=int)
args = parser.parse_args()

fpath = 'scratch'
D = np.random.randint(2,5)
lb = np.array([-1.]*D)
ub = np.array([1.]*D)
while True:
    #lengthscales from 0.05 to 1.5
    lengthscale = [int(np.random.uniform(0.05,1.4)*100)/100. for i in range(D)]
    noisevar = int(np.random.uniform(-10,-2)*10)/10.
    fnamebase = 'eihyp_d{}_n{}_l'.format(D,int(-10*noisevar)) + ''.join(['_{}'.format(int(100*i)) for i in lengthscale])
    if not any([i.startswith(fnamebase) for i in os.listdir(fpath)]):
        break
maxit = 15
#outputscale will be normalized to 1
for i in range(6):
    f, xm, truemin = objectives.genmat52ojf(D,lb,ub,A=1.,ls=lengthscale,fixs=10**noisevar,ki=GPdc.SQUEXP)
    fname = fnamebase+'_m{}_j{}.csv'.format(maxit,i)

    C=gpbo.core.config.eihypdefault(f,D,maxit,10**noisevar,fpath,fname)
    C.reccpara['kindex']=C.aqpara['kindex']=GPdc.SQUEXP
    C.reccpara['mprior']=C.aqpara['mprior']= np.array([0.]+[-1.]*D)
    C.reccpara['sprior']=C.aqpara['sprior']= np.array([0.5]*(D+1))
    C.aqpara['DH_SAMPLES']=200

    out = gpbo.search(C)