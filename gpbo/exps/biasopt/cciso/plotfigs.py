import gpbo
import numpy as np
import scipy as sp

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


mode='plot'
#mode='plot'
vers=[2,3][0]

nreps=5
D=2

s=1e-5
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objectivecciso import f

#from objective import truemin
all2confs=[]
all3confs=[]
rpath='results'

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

all2confs.append(['pesbs',C])
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
all3confs.append(['fabmod2',C])

#fabolas
C={'ninit':20,
   'nsteps':60}
#all3confs.append(['fabolas',C])


axisset={12:[1,100000,1e-4,10.],13:[1,1e5,1e-4,10.]}
labelfn = lambda x: {'eimle':'EI','pesfs':'PES','pesbs_ls':'EnvPES2','pesbs':'EnvPES','fabmod2':'FabolasM'}[x]
gpbo.plotall(all2confs+all3confs,9,rpath,logx=True,labelfn=labelfn,allylabel='Pseudo-Regret',needed=[6,13],legend=True,trueopt=10.5253080432,axisset=axisset)
gpbo.plotall(all2confs+all3confs,9,rpath,logx=True,labelfn=labelfn,allylabel='Pseudo-Regret',needed=[4,5,12],legend=False,trueopt=10.525308432,axisset=axisset)
