import gpbo
import numpy as np
import scipy as sp

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


mode=['run','plot'][1]
#mode='plot'
vers=[2,3][1]

nreps=4
D=2

s=1e-5
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

#from objective import truemin
all2confs=[]
all3confs=[]
rpath='icmlF6'

#eimle-------------------------------
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=C.reccpara['onlyafter']=10
C.stoppara = {'tmax': 60*60*10}
C.stopfn = gpbo.core.optimize.totaltstopfn

all2confs.append(['eimle',C])
#pesfs-------------------------------
C=gpbo.core.config.pesfsdefault(f,D,12,s,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=C.reccpara['onlyafter']=10

all2confs.append(['pesfs',C])


#pesfs---------------------------------
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
C.aqpara['nrandinit']=C.reccpara['onlyafter']=20


all2confs.append(['pesbs',C])

#fabolas----------------------------------
C={'ninit':20,
   'nsteps':200}
all3confs.append(['fabolas',C])

#---------------
#fabolas
C={'ninit':30,
   'nsteps':80,
   'switchkernel':True,
   'switchestimator':True}
#all3confs.append(['fabmod',C])

labelfn = lambda x: {'eimle':'EI','pesfs':'PES','pesbs':'EnvPES','fabolas':'Fabolas','fabmod':'FabolasM'}[x]
axisset={12:[30,3*1e4,2*1e-2,1],13:[0.25*1e2,4*1e4,2*1e-2,1]}
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,7,rpath,trueopt=1e-99,logx=True,labelfn=labelfn,skipinit=True,axisset=axisset,thirteenylabel='Classifier Error')
else:
    pass
