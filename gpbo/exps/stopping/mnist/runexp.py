import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

mode='run'
nreps=1
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)
parser.add_argument('-e', '--expnum', dest='expnum', action='store', default=-1, type=int)

args = parser.parse_args()

vers=2
D=2

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


from objective import f

truemin =0.
all2confs=[]
all3confs=[]
rpath='results'


#eimle-------------------------------
C=gpbo.core.config.eimledefault(f,D,10,s,rpath,'null.csv')
C.aqpara['nrandinit']=C.reccpara['onlyafter']=10
C.stoppara = {'nmax': 120}
C.stopfn = gpbo.core.optimize.nstopfn

all2confs.append(['eimle',C])

#pesfs-------------------------------
C=gpbo.core.config.pesfsdefault(f,D,10,s,rpath,'null.csv')
C.stoppara = {'nmax': 120}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['nrandinit']=C.reccpara['onlyafter']=10

all2confs.append(['pesfs',C])


#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,120,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-2
all2confs.append(['switching_2',C])

#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,120,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-4
all2confs.append(['switching_4',C])

#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,120,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-6
all2confs.append(['switching_6',C])

if args.expnum == -1:
    confs2 = all2confs
else:
    confs2 = [all2confs[args.expnum]]

gpbo.runexp(f, lb, ub, rpath, nreps, confs2, indexoffset=args.offset * nreps)
