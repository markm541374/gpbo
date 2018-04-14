import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)
parser.add_argument('-r', '--regret', dest='regret', action='store', default=-3,type=int)

args = parser.parse_args()


D=2

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)

from objective import f
def f(x,**ev):
    return -np.cos(2*x[0])-np.cos(2*x[1]),1.,{}
rawf = objectives.camel3
#def f(x,**kwargs):
#    y,c,aux = rawf(x,**kwargs)
#    return np.log(y+1),c,aux
truemin =0.
rpath='results'

#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,300,s,rpath,'switchingp_{}_{}.csv'.format(args.regret,args.offset))

C.choosepara['regretswitch']=10**args.regret

out = gpbo.search(C)
print(out)
