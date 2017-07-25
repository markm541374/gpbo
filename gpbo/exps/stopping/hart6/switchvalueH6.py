import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

mode=['run','plot'][1]
nreps=1
import argparse
gpbo.core.debugoutput['logstate']=True
gpbo.core.debugoutput['adaptive']=False
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)
parser.add_argument('-e', '--expnum', dest='expnum', action='store', default=-1, type=int)
args = parser.parse_args()

vers=[2,3][0]
D=6

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


f = objectives.shifthart6
truemin =0.
all2confs=[]
all3confs=[]
rpath='v2'
#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'null.csv')
C.chooser=gpbo.core.choosers.always0
all2confs.append(['switching_no',C])

#-----------------------
C=gpbo.core.config.switchdefault(f,D,20,61,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-2
C.aqpara[0]['nrandinit']=C.choosepara['onlyafter']=C.reccpara[0]['onlyafter']=20
all2confs.append(['switching_2Y',C])


#t = gpbo.search(C,initdata='dbout/60.p')

#-----------------------
C=gpbo.core.config.switchdefault(f,D,30,240,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-4
all2confs.append(['switching_4',C])

#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-6
all2confs.append(['switching_6',C])
if args.expnum == -1:
    confs2 = all2confs
else:
    confs2 = [all2confs[args.expnum]]
if mode == 'run':
    if vers==2:
        gpbo.runexp(f, lb, ub, rpath, nreps, confs2, indexoffset=args.offset * nreps)
    else:
        gpbo.runexp(f, lb, ub, rpath, nreps, all3confs, indexoffset=args.offset * nreps)
elif mode=='plot':
    gpbo.plotall(all2confs + all3confs, 5, rpath, trueopt=truemin + 1e-99, logx=False)
else:
    pass

