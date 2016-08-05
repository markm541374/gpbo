import optimize
import acquisitions
import reccomenders
import objectives

import scipy as sp
import os
import sys
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.DEBUG)

sys.path.append('configs')



c=1.
s=1e-9
d=[sp.NaN]

def f(x,xa):
    u=0.5*(x[0]+x[1])
    v=0.5*(0.6*x[0]-1.4*x[1])
    return -sp.cos(u*2.8*sp.pi) - 2*sp.cos(v*1.5*sp.pi)+0.1*xa*(xa-0.8)

def g(x,xa):
    if x[0]>0 and x[1]>0:
        return 2*x[0]
    elif x[0]<=0 and x[1]<=0:
        return x[1]+x[0]
    elif x[0]<=0 and x[1]>0:
        return -x[0]+x[1]
    else:
        return 0
    return 
def plotf(f,R):
    n = 100
    x = sp.linspace(-1,1,n)
    y = sp.linspace(-1,1,n)


    Z = sp.empty([n,n])
    for i in xrange(n):
	for j in xrange(n):
		Z[i,j] = f([y[j],x[i]],0.)

    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.contour(x,y,Z,40)
    for i in xrange(R.shape[0]):
        ax.plot(R[i,0],R[i,1],'ro')
    plt.show()
    return

R = sp.zeros([100,3])
for i in xrange(100):
    R[i,0] = sp.random.uniform()*2.-1
    R[i,1] = sp.random.uniform()*2.-1
nq=25

plotf(f,R[:nq,:])
O=optimize.optstate()
for i in xrange(nq):
    x = sp.array([R[i,0],R[i,1]])
    xa = sp.array(R[i,2])
    y = f(x,xa)
    O.update(x,{'s':s,'d':d,'xa':0.},y,c)
aqfn,aqpara = acquisitions.PESbs
aqpara['lb']=[-1.,-1.]
aqpara['ub']=[1.,1.]
aqpara['ev']['s']=1e-12

print aqfn(O,{'n':nq,'d':len(aqpara['ub'])},**aqpara)