import gpbo
import numpy as np
import scipy as sp
#mode='run'
mode=['run','plot'][1]
vers=[2,3][0]
D=2
nreps=4
s=1e-6
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='results3'
#-----------------------
#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'tmax': 60*60}
C.stopfn = gpbo.core.optimize.totaltstopfn

all2confs.append(['eimle_pmin',C])

#----------------------

#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'tmax': 60*60}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.reccfn=gpbo.core.reccomenders.argminrecc

all2confs.append(['eimle_best',C])

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,4,rpath,trueopt=truemin)
else:
    pass
