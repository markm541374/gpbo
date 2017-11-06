import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


D=2

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)

f = objectives.camel6

truemin =0.
rpath='results'
g = lambda x: f(x,**{})[0]
#res = sp.optimize.minimize(g,[0.0615568613,	0.0519812367,	-0.0531377095,	0.0482651495],method='CG',tol=0.5,options={'maxiter':10})

#print(res)
#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,250,s,rpath,'s_{}.csv'.format(args.offset))
C.choosepara['regretswitch']=0.1
#
out = gpbo.search(C)
#print(out)

#confs = [['switching0',gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'s_{}.csv'.format(i))] for i in range(16)]
#gpbo.plotall(confs,10,rpath,trueopt=truemin+1e-99,logx=False,showends=True)