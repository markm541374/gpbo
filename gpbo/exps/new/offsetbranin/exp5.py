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

nreps=8
D=2

s=1e-6
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='results5'



#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 2}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='none'
all2confs.append(['pesbs_none',C])


#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 2}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
all2confs.append(['pesbs_last',C])


#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,60,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 60 * 2}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
all2confs.append(['pesbs_predict',C])




if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,8,rpath,trueopt=truemin)
else:
    pass

