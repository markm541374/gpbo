import gpbo
import numpy as np
import scipy as sp
#mode='run'
gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':False,'ctaq':False,'support':False,'adaptive':False,'acqfn1d':True}
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

mode=['run','plot'][1]
vers=[2,3][0]
D=2
nreps=1
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='exp2'
#----------------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*1}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='predict'

all2confs.append(['pesvs_p1',C])

#----------------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*3}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='predict'

all2confs.append(['pesvs_p3',C])

#----------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*5}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='predict'

all2confs.append(['pesvs_p5',C])
#----------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*8}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='predict'

all2confs.append(['pesvs_p8',C])
#----------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*15}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='predict'

all2confs.append(['pesvs_p15',C])

#pesbs----------------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*15}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='last'

all2confs.append(['pesvs_l',C])

#pesbs----------------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*15}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='none'

all2confs.append(['pesvs_n',C])


labelfn = lambda x: {'pesvs':'pesvs'}[x]
axisset={11:[0,100,1e-7,1e2]}

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,4,rpath,trueopt=truemin,logx=True)
else:
    pass
