import gpbo
from gpbo.core import objectives
import numpy as np
import scipy as sp
#gpbo.core.debugoutput['tmp']=True
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--offset', dest='offset', action='store', default=0,type=int)

args = parser.parse_args()


D=4

s=0.
lb = sp.array([-1.]*D)
ub = sp.array([1.]*D)

rawf = objectives.shifthart4

truemin =0.
rpath='results'
def f(x,**kwargs):
    y,c,aux = rawf(x,**kwargs)
    return np.log(y+1),c,aux

def q(x,**kwargs):
    t  = 1.*sp.pi/6.
    z0 = x[0]*sp.cos(t)-x[1]*sp.sin(t)
    z1 = x[0]*sp.sin(t)+x[1]*sp.cos(t)
    y  = 0.1*z0**2 + 25.*z1**2
    return np.log(y+1),1.,dict()
#f([0.78547/5. -1]*3,**{})
def g(x,y):
    return rawf(np.array(x),**{})[0],0
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
C.aqpara[0]['nrandinit']=C.reccpara[0]['onlyafter']=C.choosepara['onlyafter']=80
C.choosepara['regretswitch']=1e-2
#
#C = gpbo.core.config.eihypgamma(f,D,200,s,rpath,'S_{}.csv'.format(args.offset))
out = gpbo.search(C)
#print(out)

