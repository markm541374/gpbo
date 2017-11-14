import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput['adaptive']=False
gpbo.core.debugoutput['logstate']=False
mode=['run','plot'][0]
nreps=10
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=2

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


f = objectives.shiftbraninojf
truemin =0.
all2confs=[]
all3confs=[]
rpath='results'
#-----------------------
#eimle
C=gpbo.core.config.switchdefault(f,D,10,180,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-4
#all2confs.append(['switching_6',C])


C=gpbo.core.config.directdefault(f,D,5000,s,rpath,'null.csv')
#all2confs.append(['switching_direct',C])
all2confs.append(['switching_cmaes',C])
#out = gpbo.search(C,initdata='dbout/79.p')

C=gpbo.core.config.eimledefault(f,D,100,s,rpath,'null.csv')
C.stopfn = gpbo.optimize.nstopfn
C.stoppara['nmax']=100
#all2confs.append(['eimle_fix',C])

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,8,rpath,trueopt=truemin,logx=True)
else:
    pass
