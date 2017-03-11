import gpbo
import numpy as np
import scipy as sp
#mode='run'
import copy

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

mode=['run','plot'][1]
vers=[2,3][0]
D=2
nreps=8
lb = sp.array([-1.,-1.])
ub = sp.array([1.,1.])

from objective import f

from objective import truemin
all2confs=[]
all3confs=[]
rpath='pesnoise'

#-----------------------
#eimle
def fn(x,**ev):
    y,c = f(x,**ev)
    y+=sp.random.normal()*1.
    return y,c
C=gpbo.core.config.pesfsdefault(fn,D,12,1,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 60}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle_0',C])
#-----------------------
#eimle
def fn(x,**ev):
    y,c = f(x,**ev)
    y+=sp.random.normal()*0.1
    return y,c
C=gpbo.core.config.pesfsdefault(fn,D,12,1e-2,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 60}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle_2',C])

#----------------------
def fn(x,**ev):
    y,c = f(x,**ev)
    y+=sp.random.normal()*0.01
    return y,c

C=gpbo.core.config.pesfsdefault(fn,D,12,1e-4,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 60}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle_4',C])
#-----------------------
#eimle

def fn(x,**ev):
    y,c = f(x,**ev)
    y+=sp.random.normal()*0.001
    return y,c

C=gpbo.core.config.pesfsdefault(fn,D,12,1e-6,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 60}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle_6',C])

#-----------------------
#eimle
def fn(x,**ev):
    y,c = f(x,**ev)
    y+=sp.random.normal()*0.0001
    return y,c
C=gpbo.core.config.pesfsdefault(fn,D,12,1e-8,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 60}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle_8',C])
#-----------------------
#eimle
def fn(x,**ev):
    y,c = f(x,**ev)
    y+=sp.random.normal()*0.00001
    return y,c
C=gpbo.core.config.pesfsdefault(fn,D,12,1e-10,rpath,'null.csv')
C.aqpara['nrandinit']=10
C.stoppara = {'nmax': 60}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle_10',C])
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,6,rpath,trueopt=truemin)
else:
    pass
