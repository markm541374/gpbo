import gpbo
import numpy as np
import scipy as sp

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


mode=['run','plot'][0]
#mode='plot'
vers=[2,3][0]

nreps=1
D=2

s=1e-6
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

#from objective import truemin
all2confs=[]
all3confs=[]
rpath='results0'

#eimle-------------------------------
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 80}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle',C])

#pesbs---------------------------------
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'nmax': 140}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=20


all2confs.append(['pesbs',C])

#fabolas----------------------------------
C={'ninit':20,
   'nsteps':140}
all3confs.append(['fabolas',C])



if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,8,rpath)
else:
    pass
