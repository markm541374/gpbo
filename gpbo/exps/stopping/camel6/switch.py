import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

#gpbo.core.debugoutput=True
#gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':False,'ctaq':False,'support':False,'adaptive':True,'logstate':False}
mode='plot'
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


f = objectives.shifthart3
truemin =0
all2confs=[]
all3confs=[]
rpath='results'
#-----------------------

#all2confs.append(['switching_no',None])

all2confs.append(['switchingp_-2',None])

all2confs.append(['switchingp_-4',None])

all2confs.append(['pesfs',None])
all2confs.append(['eihyp',None])
if mode=='run':
    if vers==2:
        gpbo.runexp(f,lb,ub,rpath,nreps,all2confs,indexoffset=args.offset*nreps)
    else:
        gpbo.runexp(f,lb,ub,rpath,nreps,all3confs,indexoffset=args.offset*nreps)
elif mode=='plot':
    gpbo.plotall(all2confs+all3confs,16,rpath,trueopt=truemin+1e-99,logx=False,showends=True)
else:
    pass
