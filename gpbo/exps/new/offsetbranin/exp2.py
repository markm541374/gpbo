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

nreps=1
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
C.stoppara = {'nmax': 50 }
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['overhead']='none'
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_LAPAPROT]
C.aqpara['DM_SLICELCBPARA']=20
#all2confs.append(['pesbs_lap',C])

#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'nmax': 50 }
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['overhead']='none'
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_SLICEEI]
C.aqpara['DM_SLICELCBPARA']=2.

all2confs.append(['pesbs_ei',C])

#------------------
#pesfs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'nmax': 50 }
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['overhead']='none'
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_SLICELCB]
C.aqpara['DM_SLICELCBPARA']=2.
#all2confs.append(['pesbs_lcb',C])

#------------------
#pesfs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'nmax': 50 }
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['overhead']='none'
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_SLICEPM]
C.aqpara['DM_SLICELCBPARA']=2.
#all2confs.append(['pesfs_pm',C])


if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,1,rpath,trueopt=truemin)
else:
    pass

