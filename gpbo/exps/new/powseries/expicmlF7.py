import gpbo
import numpy as np
import scipy as sp

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


mode=['run','plot'][1]
#mode='plot'
vers=[2,3][0]

nreps=2
D=2

s=1e-5
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

#from objective import truemin
all2confs=[]
all3confs=[]
rpath='icmlF7'

#-----------------------
#eimle
C=gpbo.core.config.eimledefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=C.reccpara['onlyafter']=10
C.stoppara = {'tmax': 60*60*8}
C.stopfn = gpbo.core.optimize.totaltstopfn

all2confs.append(['eimle',C])

#pesfs
C=gpbo.core.config.pesfsdefault(f,D,12,s,rpath,'null.csv')
C.aqpara['nrandinit']=C.reccpara['onlyafter']=10
C.stoppara = {'tmax': 60*60*8}
C.stopfn = gpbo.core.optimize.totaltstopfn

all2confs.append(['pesfs',C])

#pesbs----------------------------
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 8}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['traincfn']='llogfull'
C.aqpara['overhead']='predict'
C.aqpara['hyp_chains']=6
C.aqpara['nrandinit']=C.reccpara['onlyafter']=20
C.aqpara['mprior']=C.reccpara['mprior']= sp.array([1.]+[0.]*(D+1)+[-3])
C.aqpara['sprior']=C.reccpara['sprior']= sp.array([1.]*(D+2)+[3])
C.aqpara['kindex']=C.reccpara['kindex']=gpbo.core.GPdc.MAT52CS

all2confs.append(['pesbs_ls',C])


#-----------------
#mtbo
C={'lowtask':4,
   'ninit':15,
   'nsteps':40}

#all3confs.append(['mtbo4',C])

#-----------------
#mtbo
C={'lowtask':16,
   'ninit':15,
   'nsteps':40}

#all3confs.append(['mtbo16',C])

#-----------------
#mtbo
C={'lowtask':64,
   'ninit':15,
   'nsteps':80}

#all3confs.append(['mtbo64',C])
#---------------
#fabolas
C={'ninit':20,
   'nsteps':40}
all3confs.append(['fabmod',C])

#fabolas
C={'ninit':20,
   'nsteps':60}
#all3confs.append(['fabolas',C])


labelfn = lambda x: {'eimle':'EI','pesfs':'PES','pesbs_ls':'EnvPES','fabmod':'FabolasM'}[x]
axisset={6:[1e3,4*1e4,10.5,14]}
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,8,rpath,logx=True,labelfn=labelfn,axisset=axisset,sixylabel='GP Negative Log-Likelihood')
else:
    pass
