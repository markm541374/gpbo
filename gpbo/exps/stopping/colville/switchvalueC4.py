import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#mode='run'

gpbo.core.debugoutput['adaptive']=False
gpbo.core.debugoutput['logstate']=False
gpbo.core.debugoutput['forceNoiseFloor']=-9
mode=['run','plot'][0]
nreps=10
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()

vers=[2,3][0]
D=4

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)


f = objectives.colville


truemin =0.
all2confs=[]
all3confs=[]
rpath='results'


#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,200,s,rpath,'s2.csv')
C.choosepara['regretswitch']=1e-2
#all2confs.append(['switching_6',C])

out = gpbo.search(C)
print(out)

