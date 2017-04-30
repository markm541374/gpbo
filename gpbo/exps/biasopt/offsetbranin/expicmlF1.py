import gpbo
import numpy as np
import scipy as sp
#mode='run'

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

mode=['run','plot'][1]
vers=[2,3][0]
D=2
nreps=1
s=1e-6
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='F1new'
#pesbs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'nmax': 100}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['nrandinit']=10
C.reccfn=gpbo.core.reccomenders.argminrecc

all2confs.append(['pesfs_argmin',C])
#----------------------#pesbs
C=gpbo.core.config.pesfsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'nmax': 100}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['nrandinit']=10


all2confs.append(['pesfs_predmin',C])

#-----------------------
#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 100}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle_predmin',C])

#----------------------

#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 100}
C.stopfn = gpbo.core.optimize.nstopfn
C.reccfn=gpbo.core.reccomenders.argminrecc

all2confs.append(['eimle_argmin',C])
labelfn = lambda x: {'eimle_argmin':'EI-AM','pesfs_predmin':'PES-PM','pesfs_argmin':'PES-AM','eimle_predmin':'EI-PM','fabmod':'FabolasM'}[x]
axisset={11:[0,150,1e-8,1e2]}

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,25,rpath,trueopt=truemin,labelfn=labelfn,axisset=axisset)
else:
    pass
