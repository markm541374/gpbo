import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput['adaptive']=False
gpbo.core.debugoutput['logstate']=False
mode=['run','plot'][0]
nreps=25
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=3

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


f = objectives.shifthart3


truemin =0.
all2confs=[]
all3confs=[]
rpath='eipi'


#-----------------------
C=gpbo.core.config.eihypdefault(f,D,160,s,rpath,'null.csv')
C.stopfn = gpbo.optimize.PIorNstopfn
C.stoppara['PImin']=0.1
#all2confs.append(['eihyp_01',C])


#-----------------------
C=gpbo.core.config.eihypdefault(f,D,160,s,rpath,'null.csv')
C.stopfn = gpbo.optimize.PIorNstopfn
C.stoppara['PImin']=0.01
#all2confs.append(['eihyp_001',C])

#-----------------------
C=gpbo.core.config.eihypdefault(f,D,160,s,rpath,'null.csv')
C.stopfn = gpbo.optimize.PIorNstopfn
C.stoppara['PImin']=0.001
#all2confs.append(['eihyp_0001',C])

#-----------------------
C=gpbo.core.config.pesfsdefault(f,D,160,s,rpath,'null.csv')
C.stopfn = gpbo.optimize.PIorNstopfn
C.stoppara['PImin']=0.1
all2confs.append(['pesfs_01',C])


#-----------------------
C=gpbo.core.config.pesfsdefault(f,D,160,s,rpath,'null.csv')
C.stopfn = gpbo.optimize.PIorNstopfn
C.stoppara['PImin']=0.01
all2confs.append(['pesfs_001',C])

#-----------------------
C=gpbo.core.config.pesfsdefault(f,D,160,s,rpath,'null.csv')
C.stopfn = gpbo.optimize.PIorNstopfn
C.stoppara['PImin']=0.001
all2confs.append(['pesfs_0001',C])
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,25,rpath,trueopt=truemin+1e-99,logx=False,showends=True)
else:
    pass
