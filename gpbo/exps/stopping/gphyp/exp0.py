import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':False,'ctaq':False,'support':False,'adaptive':True,'logstate':False}
mode='plot'
nreps=1
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=3

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


f = objectives.shifthart3
truemin =0.0145
all2confs=[]
all3confs=[]
rpath='r2'
#-----------------------

#all2confs.append(['switching_no',None])

all2confs.append(['switching_2',None])

all2confs.append(['switching_4',None])

all2confs.append(['switching_6',None])

all2confs.append(['eimle',None])

all2confs.append(['pesfs',None])

gpbo.plotall(all2confs+all3confs,3,rpath,trueopt=truemin+1e-99,logx=False,showends=True)
