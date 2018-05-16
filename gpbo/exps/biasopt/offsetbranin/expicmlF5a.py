import gpbo
import numpy as np
import scipy as sp
#mode='run'

mode=['run','plot'][1]
nreps=1
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=2

s=1e-3
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='F5anew'
#-----------------------
#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'tmax': 60*60*15}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.reccfn = gpbo.core.reccomenders.argminrecc
all2confs.append(['eimle',C])

#pesfs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 15}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=10
all2confs.append(['pesfs',C])

#pesbs----------------------------
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 15}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
C.aqpara['nrandinit']=20

all2confs.append(['pesbs',C])

#-----------------
#mtbo
C={'lowtask':4,
   'ninit':20,
   'switchestimator':True,
   'nsteps':100}

#all3confs.append(['mtbo4',C])

#-----------------
#mtbo
C={'lowtask':64,
   'ninit':20,
   'switchestimator':True,
   'nsteps':150}

#all3confs.append(['mtbo64',C])
#---------------
#fabolas
C={'ninit':30,
   'nsteps':80}
#all3confs.append(['fabolas',C])

#---------------
#fabolas
C={'ninit':20,
   'nsteps':80,
   'switchkernel':True,
   'switchestimator':True}
all3confs.append(['fabmod',C])

#--------------

labelfn = lambda x: {'eimle':'EI','pesfs':'PES','pesbs':'EnvPES','fabolas':'Fabolas','fabmod':'FabolasM'}[x]
#axisset={12:[1000,60*60*20,1e-8,100],13:[1000,60*60*20,1e-8,100]}
axisset={12:[1e3,1e5,1e-8,100],13:[1e3,1e5,1e-8,100]}
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,20,rpath,trueopt=truemin,labelfn=labelfn,logx=True,axisset=axisset,needed=[12],legend=True)
    gpbo.plotall(all2confs+all3confs,20,rpath,trueopt=truemin,labelfn=labelfn,logx=True,axisset=axisset,needed=[13],legend=False)
else:
    pass
