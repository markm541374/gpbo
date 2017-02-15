import gpbo
import numpy as np
import scipy as sp
#mode='run'

mode=['run','plot'][1]
nreps=8
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=3

s=1e-8
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='results1'
#-----------------------
#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'tmax': 60*60*20}
C.stopfn = gpbo.core.optimize.totaltstopfn
all2confs.append(['eimle',C])


#pesfs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 20}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
C.aqpara['nrandinit']=10
all2confs.append(['pesfs',C])

#pesbs----------------------------
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 20}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
C.aqpara['nrandinit']=20

all2confs.append(['pesbs',C])

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
   'nsteps':200,
   'switchkernel':True,
   'switchestimator':True}
#all3confs.append(['fabmod',C])
#---------------
#fabolas
C={'ninit':20,
   'nsteps':140}
#all3confs.append(['fabolas',C])
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,8,rpath,trueopt=truemin)
else:
    pass
