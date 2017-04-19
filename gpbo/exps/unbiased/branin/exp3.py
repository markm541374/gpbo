import gpbo
import numpy as np
import scipy as sp
#mode='run'
gpbo.core.debugoutput=True
from collections import defaultdict
gpbo.core.debugoptions=defaultdict(lambda :False)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

mode=['run','plot'][0]
vers=[2,3][0]
D=2
nreps=8
s=1e-12
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='exp3'

#---------------------------------------
C=gpbo.core.config.eimledefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'nmax': 200}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['nrandinit']=10
C.aqpara['smode']=C.reccpara['smode']='direct'
all2confs.append(['eimle_direct',C])

#---------------------------------------
C=gpbo.core.config.eimledefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'nmax': 200}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['nrandinit']=10
C.aqpara['smode']=C.reccpara['smode']='multi'
all2confs.append(['eimle_multi',C])

#---------------------------------------
C=gpbo.core.config.eimledefault(f,D,50,s,rpath,'null.csv')
C.stoppara = {'nmax': 200}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['nrandinit']=10
C.aqpara['smode']=C.reccpara['smode']='dthenl'
all2confs.append(['eimle_both',C])

labelfn = lambda x: {'pesvs':'pesvs'}[x]
axisset={11:[0,100,1e-7,1e2]}

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,8,rpath,trueopt=truemin,logx=True)
else:
    pass
