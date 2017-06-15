import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput['adaptive']=False
gpbo.core.debugoutput['logstate']=False
mode=['run','plot'][0]
nreps=1
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=3

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)



f = objectives.shifthart3D
truemin =0.
all2confs=[]
all3confs=[]
rpath='results'



#-----------------------
C=gpbo.core.config.eihypdefault(f,D,10,s,'results','null.csv')
C.stoppara = {'nmax': 40}
C.stopfn = gpbo.core.optimize.nstopfn
C.ojfchar['batchgrad']=True
#all2confs.append(['eihypgrad',C])

#-----------------------
C=gpbo.core.config.pesfsdefault(f,D,10,s,'results','null.csv')
C.stoppara = {'nmax': 40}
C.stopfn = gpbo.core.optimize.nstopfn
C.ojfchar['batchgrad']=True
all2confs.append(['pesfsgradx',C])

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,14,rpath,trueopt=truemin+1e-99,logx=False,showends=True)
else:
    pass
