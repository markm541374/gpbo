import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
gpbo.core.debugoutput['tmp']=True
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


D=2

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)

rawf = objectives.camel6

truemin =0.
rpath='results'
def f(x,**kwargs):
    y,c,aux = rawf(x,**kwargs)
    return np.log(y+1),c,aux

g = lambda x: f(x,**{})[0]
#from matplotlib import pyplot as plt
#m=200
#x = np.linspace(-1,1,m)
#y = np.linspace(-1,1,m)
#A = np.empty([m,m])
#for i in range(m):
#    for j in range(m):
#        A[j,i] = np.log(g([x[i],y[j]])+1)
#CS = plt.contour(x,y,A,40)

#res = sp.optimize.minimize(g,[0.0615568613,	0.0519812367,	-0.0531377095,	0.0482651495],method='CG',tol=0.5,options={'maxiter':10})

#print(res)
#-----------------------
C=gpbo.core.config.switchdefault(f,D,10,250,s,rpath,'s_{}.csv'.format(args.offset))
C.choosepara['regretswitch']=1e-4
#
out = gpbo.search(C)
#print(out)

#confs = [['switching0',gpbo.core.config.switchdefault(f,D,10,160,s,rpath,'s_{}.csv'.format(i))] for i in range(16)]
#gpbo.plotall(confs,10,rpath,trueopt=truemin+1e-99,logx=False,showends=True)