import gpbo
import numpy as np
import scipy as sp
#mode='run'
mode=['run','plot'][1]
vers=[2,3][0]
D=2
nreps=10
s=1e-6
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='results3x'
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

all2confs.append(['eimle_pmin',C])

#----------------------

#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 100}
C.stopfn = gpbo.core.optimize.nstopfn
C.reccfn=gpbo.core.reccomenders.argminrecc

all2confs.append(['eimle_best',C])
labelfn = lambda x: {'eimle_best':'EI-AM','pesfs_predmin':'PES_PM','pesfs_argmin':'PES_AM','eimle_pmin':'EI-PM','fabmod':'FabolasM'}[x]

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,10,rpath,trueopt=truemin,labelfn=labelfn)
else:
    pass
