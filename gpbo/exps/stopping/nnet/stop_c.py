import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
import nnetfns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)
parser.add_argument('-r', '--regret', dest='regret', action='store', default=-2,type=int)

args = parser.parse_args()


D=4

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)

#rawf = objectives.shifthart4
#def f(x,**kwargs):
#    y,c,aux = rawf(x,**kwargs)
#    return np.log(y+1),c,aux
f=nnetfns.nnet_cancer
#truemin =0.
rpath='results'

#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,300,s,rpath,'switchingp_{}_{}.csv'.format(args.regret,args.offset))

C.choosepara['regretswitch']=10**args.regret

out = gpbo.search(C)
print(out)
