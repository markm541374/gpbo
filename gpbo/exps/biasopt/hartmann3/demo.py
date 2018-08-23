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

s=1e-8
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)

from objective import f

#from objective import truemin
truemin = -3.32237
all2confs=[]
all3confs=[]
rpath='tmp'
#-----------------------
#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'tmax': 60*60*10}
C.stopfn = gpbo.core.optimize.totaltstopfn
#all2confs.append(['eihyp',C])


#pesfs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
C.aqpara['nrandinit']=10
all2confs.append(['pesfs_-4',C])

#pesbs----------------------------
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
C.aqpara['nrandinit']=20

#all2confs.append(['pesbs',C])

#-----------------
#mtbo
C={'lowtask':4,
   'ninit':20,
   'nsteps':80}

#all3confs.append(['mtbo4',C])

#-----------------
#mtbo
C={'lowtask':16,
   'ninit':20,
   'nsteps':200,
   'switchestimator':True}

#all3confs.append(['mtbo16',C])

#-----------------
#mtbo
C={'lowtask':64,
   'ninit':20,
   'nsteps':150}

#all3confs.append(['mtbo8',C])
#---------------
#fabolas
C={'ninit':20,
   'nsteps':100,
   'switchkernel':True,
   'switchestimator':True}
#all3confs.append(['fabmod',C])
#---------------
#fabolas
C={'ninit':20,
   'nsteps':140}
#all3confs.append(['fabolas',C])
labelfn = lambda x: {'eihyp':'EI','pesfs_-4':'PES','pesbs':'EnvPES','fabmod':'FabolasM'}[x]
#axisset={12:[1000,70000,1e-6,10],13:[2000,24*60*60,1e-6,10]}
axisset={12:[1e3,1e5,1e-6,10],13:[1e3,1e5,1e-6,10]}
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,1,rpath,trueopt=0.*truemin,logx=False,labelfn=labelfn,axisset=axisset,needed=[4],legend=False,offs=9)
else:
    pass
