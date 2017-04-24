import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':False,'ctaq':False,'support':False,'adaptive':True,'logstate':True}
mode=['run','plot'][0]
nreps=1
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=2

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


f, xm, truemin = objectives.genmat52ojf(D,lb,ub,ls=0.2,fixs=-1)
all2confs=[]
all3confs=[]
rpath='2d'
#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'null.csv')
C.chooser=gpbo.core.choosers.always0
#all2confs.append(['switching_no',C])

#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-4
#C.stopfn=gpbo.core.nstopfn
#C.stoppara['n']=35
#out = gpbo.search(C,initdata='dbout/34.p')
all2confs.append(['switching_4',C])


#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-3
#all2confs.append(['switching_3',C])

#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'null.csv')
C.choosepara['regretswitch']=1e-2
#all2confs.append(['switching_2',C])


if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,1,rpath,trueopt=truemin+1e-99,logx=False)
else:
    pass
