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

mode=['run','plot'][0]
vers=[2,3][0]
D=2
nreps=8
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='exp1'
#pesbs----------------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*6}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='predict'

#all2confs.append(['pesvs_p',C])
#----------------------#pesvs
#pesbs----------------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*6}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='last'

#all2confs.append(['pesvs_l',C])

#pesbs----------------------------
C=gpbo.core.config.pesvsdefault(f,D,50,-6,0,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*6}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='none'

#all2confs.append(['pesvs_n',C])

#pesbs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,1e-7,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*6}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='none'

#all2confs.append(['pesfs_7',C])
#pesbs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,1e-5,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*6}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='none'

#all2confs.append(['pesfs_5',C])
#pesbs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,1e-3,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*6}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='none'

#all2confs.append(['pesfs_3',C])
#pesbs----------------------------
C=gpbo.core.config.pesfsdefault(f,D,50,1e-1,rpath,'null.csv')
C.stoppara = {'tmax': 60*60*6}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['nrandinit']=10
C.aqpara['overhead']='none'

#all2confs.append(['pesfs_1',C])

labelfn = lambda x: {'pesvs':'pesvs'}[x]
axisset={11:[0,100,1e-7,1e2]}

if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,1,rpath,trueopt=truemin,labelfn=labelfn,axisset=axisset)
else:
    pass
