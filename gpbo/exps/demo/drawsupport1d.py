import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from gpbo.core import objectives
import gpbo
ki = gpbo.core.GPdc.SQUEXP
lb = np.array([-1.])
ub = np.array([1.])
D=1

[X,Y,S,D] = gpbo.core.ESutils.gen_dataset(8, D, lb, ub, ki, sp.array([1,0.3]),s=1e-4)

G = gpbo.core.PES.makeG(X,Y,S,D,ki, np.array([1.,0.]),np.array([1.,1.]),36)

f,ax = plt.subplots(nrows=2, ncols=1,figsize=[10,8])

n = 160
z = np.linspace(-1,1,n)
m,v = G.infer_diag_post(z,[[np.NaN]]*n)
ax[0].plot(z,m.flatten(),'b')
ax[0].fill_between(z,m.flatten()-2*np.sqrt(v).flatten(),m.flatten()+2*np.sqrt(v).flatten(),facecolor='b',edgecolor='Aliceblue',alpha=0.2)
ax[0].plot(X,Y,'ro')

m = 250
r = np.linspace(-1,1,250)
c = np.zeros(m)
for i in range(10):
    P = G.draw_post(r,[[np.NaN]]*m,400)
    A = np.argmin(P,axis=1)
    for a in A:
        c[a]+=1
ax[1].plot(r.flatten(),c.flatten()/np.sum(c),'r')


Z = gpbo.core.ESutils.draw_support(G, lb, ub, 4000,gpbo.core.ESutils.SUPPORT_LAPAPROT,para=8,pad_unif=False)
f2,ax2 = plt.subplots()
n,bins,patches = ax2.hist(Z,bins=np.linspace(-1,1,m+1),normed=1)
ax[1].plot(bins[1:]-0.5*(bins[1]-bins[0]),n/np.sum(n),'g')


Z = gpbo.core.ESutils.draw_support(G, lb, ub, 4000,gpbo.core.ESutils.SUPPORT_SLICEEI)
n,bins,patches = ax2.hist(Z,bins=np.linspace(-1,1,m+1),normed=1)
ax[1].plot(bins[1:]-0.5*(bins[1]-bins[0]),n/np.sum(n),'b')

Z = gpbo.core.ESutils.draw_support(G, lb, ub, 4000,gpbo.core.ESutils.SUPPORT_SLICELCB,para=2)
n,bins,patches = ax2.hist(Z,bins=np.linspace(-1,1,m+1),normed=1)
ax[1].plot(bins[1:]-0.5*(bins[1]-bins[0]),n/np.sum(n),'purple')

f.savefig('figs/1ddemo.png')
