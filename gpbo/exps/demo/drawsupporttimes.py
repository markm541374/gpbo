import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from gpbo.core import objectives
import gpbo
import time
ki = gpbo.core.GPdc.SQUEXP
lb = np.array([-1.,-1.,-1.])
ub = np.array([1.,1.,1.])
Dim=3

X = np.random.uniform(-1,1,size=[50,Dim])
Y = np.empty([X.shape[0],1])
for i in range(X.shape[0]):
    Y[i,0] = gpbo.core.objectives.shifthart3(X[i,:],**{'s':1e-6,'d':[[sp.NaN]]})[0]
S = 1e-6 * np.ones([X.shape[0],1])
D = [[np.NaN]]*X.shape[0]

G = gpbo.core.PES.makeG(X,Y,S,D,ki, np.array([1.,0.,0.,0.]),np.array([1.,1.,1.,1.]),36)



def times(f,nmax):
    nt=25

    times = np.zeros([2,nt])
    for j,n in enumerate([int(i) for i in np.linspace(1,nmax,nt)]):
        print(j)
        t0=time.clock()
        try:
            Z = f(n)
        except:
            pass
        t1 = time.clock()
        times[1,j] = t1-t0
        times[0,j] = n
    return times



ei = lambda n: gpbo.core.ESutils.draw_support(G, lb, ub, n,gpbo.core.ESutils.SUPPORT_SLICEEI)
eitimes = times(ei,1000)

lcb = lambda n: gpbo.core.ESutils.draw_support(G, lb, ub, n,gpbo.core.ESutils.SUPPORT_SLICELCB,para=2)
#lcbtimes = times(ei,200)

laprot = lambda n: gpbo.core.ESutils.draw_support(G, lb, ub, n,gpbo.core.ESutils.SUPPORT_LAPAPROT,para=20,pad_unif=False)
laptimes = times(laprot,1000)

laprotw2 = lambda n: gpbo.core.ESutils.draw_support(G, lb, ub, n,gpbo.core.ESutils.SUPPORT_LAPAPROT,para=20,pad_unif=False,weighted=2)
lapw2times = times(laprotw2,1000)

import emcee
def emceeei(n):
    print('draw under emceeei')
    def f(x,G):
        if all(x>lb) and all(x<ub):
            try:
                ei=G.infer_EI(sp.array(x),[[sp.NaN]])[0,0]
                if np.isnan(np.log(ei)):
                    return 1e-99
                #print ei
                return np.log(ei)
            except:
                #ei=g.infer_EI(sp.array(x),[[sp.NaN]])[0,0]
                G.printc()
                raise
        else:
            return -1e99

    nwalkers = 24
    p0 = np.random.uniform(-1,1,Dim * nwalkers).reshape((nwalkers, Dim))
    sampler = emcee.EnsembleSampler(nwalkers, Dim, f, args=[G])
    sampler.reset()
    sampler.run_mcmc(p0, max(1,n//nwalkers))
    return sampler.flatchain
#emeitimes = times(emceeei,5000)

f,ax = plt.subplots(nrows=1, ncols=1,figsize=[8,8])
ax.plot(laptimes[0,:],laptimes[1,:],'bo')
ax.plot(lapw2times[0,:],lapw2times[1,:],'ro')
ax.plot(eitimes[0,:],eitimes[1,:],'go')
ax.set_yscale('log')
#ax.plot(lcbtimes[0,:],lcbtimes[1,:],'g-.')
#ax.plot(emeitimes[0,:],emeitimes[1,:],'g')


f.savefig('figs/timedemo.png')
