import numpy as np
import scipy as sp
import gpbo
from gpbo.core import GPdc
from gpbo.core import objectives
import argparse
import pickle
import time
from scipy.interpolate import interp2d
from scipy.optimize import minimize
from scipy.stats import gamma as gammadist
from scipy.stats import norm
from scipy.stats import lognorm
from scipy import stats
from scipy.special import gamma

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index', dest='index', action='store', default=0,type=int)
parser.add_argument('-d', '--dimension', dest='dimension', action='store', default=2,type=int)
args = parser.parse_args()

D = 2#args.dimension
lb = np.array([-1.]*D)
ub = np.array([1.]*D)
#lengthscales from 0.05 to 1.5
lengthscale = [lognorm.rvs(s=0.5*np.log(10),scale=10**-0.5) for i in range(D)]#,0.8,0.7]
#lengthscale = [0.2]*D
#outputscale will be normalized to 1
fc, xm, truemin = objectives.genmat52ojf(D,lb,ub,A=1.,ls=lengthscale,fixs=-1,ki=GPdc.SQUEXP)

def cfn(sigma):
    return 30*(1e-6)/sigma
def icfn(c):
    return 30*(1e-6)/c

def f(x,**ev):
    y,c,aux = fc(x,**ev)
    return y,cfn(ev['s']),aux

fpath = 'results'
from gpbo.core import debugoutput


vtarget = 2.3579603719055325e-06

for i in range(1):
    fname = 'eihypV_{}_{}_{}.csv'.format(int(lengthscale[0]*100),int(lengthscale[1]*100),i)
    C=gpbo.core.config.eihypdefault(f,D,10,vtarget,fpath,fname)
    C.reccpara['kindex']=C.aqpara['kindex']= GPdc.SQUEXP
    C.reccpara['mprior']=C.aqpara['mprior']= sp.array([0.]+[-0.5]*D)
    C.reccpara['sprior']=C.aqpara['sprior']= sp.array([0.5]*(D+1))

    C.aqpara['DH_SAMPLES']=200
    C.aqpara['B']=75*cfn(1e-6)
    C.stopfn = gpbo.optimize.totalTorNstopfn
    C.stoppara['nmax']=200
    C.stoppara['tmax']=C.aqpara['B']
    C.aqpara['costfn']=cfn
    C.aqpara['icostfn']=icfn


    #mean-var of regret model

    out = gpbo.search(C)
