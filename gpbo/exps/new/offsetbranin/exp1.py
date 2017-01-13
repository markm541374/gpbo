import gpbo
import numpy as np
import scipy as sp
mode=['run','plot'][1]
#mode='plot'

D=2

s=1e-6
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
allconfs=[]
rpath='results1'



#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 55}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='none'
C.aqpara['nrandinit']=16

allconfs.append(['pesbs_noov',C])

#-----------------
#pesbs
C=gpbo.core.config.pesbsdefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'tmax': 60 * 55}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='last'
C.aqpara['nrandinit']=16

allconfs.append(['pesbs_ov',C])

if mode=='run':
    gpbo.runexp(f,lb,ub,rpath,1,allconfs)
elif mode=='plot':
    gpbo.plotall(allconfs,1,rpath,trueopt=truemin)
else:
    pass
