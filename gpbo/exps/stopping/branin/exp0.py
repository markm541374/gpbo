import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput['adaptive']=False
gpbo.core.debugoutput['logstate']=True
mode=['run','plot'][0]
nreps=1
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
rpath='tmp'
#-----------------------
#eimle
C=gpbo.core.config.switchdefault(f,D,10,180,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-4
all2confs.append(['switching_6',C])

#out = gpbo.search(C,initdata='dbout/79.p')
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,8,rpath,trueopt=truemin,logx=True)
else:
    pass
