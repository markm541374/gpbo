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

nreps=2
D=2

s=1e-6
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='results2'



#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=16
C.aqpara['DH_SAMPLES']=8
C.aqpara['DM_SAMPLES']=16
C.aqpara['DM_SUPPORT']=1000
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_LAPAPROT]
C.aqpara['DM_SLICELCBPARA']=20
all2confs.append(['pesbs_v0',C])

#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=16
C.aqpara['DH_SAMPLES']=16
C.aqpara['DM_SAMPLES']=32
C.aqpara['DM_SUPPORT']=1000
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_LAPAPROT]
C.aqpara['DM_SLICELCBPARA']=20

all2confs.append(['pesbs_v1',C])

#------------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=16
C.aqpara['DH_SAMPLES']=32
C.aqpara['DM_SAMPLES']=64
C.aqpara['DM_SUPPORT']=1000
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_LAPAPROT]
C.aqpara['DM_SLICELCBPARA']=20
all2confs.append(['pesbs_v2',C])

#------------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 10}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=16
C.aqpara['DH_SAMPLES']=16
C.aqpara['DM_SAMPLES']=32
C.aqpara['DM_SUPPORT']=500
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_LAPAPROT]
C.aqpara['DM_SLICELCBPARA']=20
all2confs.append(['pesbs_v3',C])


if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,2,rpath,trueopt=truemin)
else:
    pass

