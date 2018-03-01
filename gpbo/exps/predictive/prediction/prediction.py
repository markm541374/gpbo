import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import dill as pickle
import sys
import tqdm
import os
import time
import gpbo
cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

_,OP = pickle.load(open('../overhead/overmodel.p'))
from gpbo.exps.predictive.overhead.overheads import cvmodel
def OM(X):
    return cvmodel(X, OP)

_,PP = pickle.load(open('../performance/model.p'))
from gpbo.exps.predictive.performance.model import perfmodel
def PM(X,S,L):
    X = X.reshape(-1,1)
    truemean,log10std = perfmodel(np.hstack([X,np.log10(S)*np.ones_like(X)]),L,PP)
    natlogmu = np.log(truemean)
    natlogvar = (np.log(10)*log10std)**2
    return natlogmu,natlogvar


def probsteps(B,C,ax=None):
    mx=10
    while True:
        #availabel step numbers are up to B/C
        nsteps = np.arange(max(1,min(mx,1+(B-10*C)/C)))
        if len(nsteps)==0:
            return np.arange(0),np.zeros(0)
        #overhead budget remaining at step
        bover = B-C*(nsteps+10)
        #perstep mean and cov for overhead
        pm,pv = OM(nsteps)
        #cumulative mean and std
        cmean = np.cumsum(pm)
        cstd = np.sqrt(np.array([np.sum(pv[:i+1,:i+1]) for i in range(pv.shape[0])]))
        cnstep = sp.stats.norm.cdf(bover,loc=cmean,scale=cstd)
        cnstep[np.isnan(cnstep)]=1
        cnstep = np.minimum.accumulate(cnstep)
        if cnstep[-1]<1e-6 or len(nsteps)<mx:
            break
        mx*=5

    pnstep = np.nan_to_num(np.hstack([-np.diff(cnstep),0]))
    if len(pnstep)==1:
        pnstep=np.ones(1)
    pnstep = pnstep/np.sum(pnstep)
    #initialization offset
    nsteps+=10
    if not ax is None:
        ax.plot(nsteps,cmean,cols[0])
        ax.plot(nsteps,cmean-2*cstd,cols[0],linestyle='--')
        ax.plot(nsteps,cmean+2*cstd,cols[0],linestyle='--')

        ax.plot([0.,B/C],[B,0.],cols[0],linewidth=0.5)
        axt = ax.twinx()
        axt.plot(nsteps,cnstep,cols[1],linewidth=0.75)
    #if np.any(np.isnan(pnstep)):
    #    pass
    return nsteps, pnstep

def muss2mv(mu,ss):
    #natural log
    m = np.exp(mu+0.5*ss)
    v = (np.exp(ss)-1)*np.exp(2*mu+ss)
    return m,v
def mv2muss(m,v):
    #natural log
    mu = np.log(m/np.sqrt(1+v/m**2))
    ss = np.log(1+v/m**2)
    return mu,ss

def perfatl(nsteps,probn,S,L,ax=None):
    #probn = np.zeros_like(probn)
    #probn[10]=1.
    mu,ss = PM(nsteps,S,L)
    #lognormal mean variance in natural log

    m,v = muss2mv(mu,ss)

    margm = np.sum(probn*m)
    margv = np.sum(probn*(m-margm)**2)+np.sum(probn*v)
    margmu,margss = mv2muss(margm,margv)


    if not ax is None:
        ax.plot(nsteps,np.exp(mu),cols[0])
        ax.plot(nsteps,np.exp(mu+2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(nsteps,np.exp(mu-2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(nsteps,np.exp(mu+0.5*ss),cols[0],linestyle='-.')

        axt = ax.twinx()
        axt.plot(nsteps,probn,cols[1],linewidth=0.75)

        ax.set_yscale('log')
        iplt = np.argmax(probn)
        xplt = nsteps[iplt]
        ax.plot(xplt,margm,'rx')
        ax.plot(xplt,np.exp(margmu),'r.')
        ax.plot([xplt,xplt],[np.exp(margmu-2*np.sqrt(margss)),np.exp(margmu+2*np.sqrt(margss))],'r')
    return margmu, margss

def steps(nsteps,probn):
    #returns mean and 95% center interval step numbers
    mean = np.sum(probn*nsteps)
    if len(nsteps)==0:
        return 10,10,10
    lower = nsteps[np.argmax(np.cumsum(probn)>0.025)]
    upper = nsteps[np.argmax(np.cumsum(probn)>0.975)]
    return lower,mean,upper

def perfatBCVL(B,C,V,L):
    nsteps,probn = probsteps(B,C)
    mu,ss = perfatl(nsteps,probn,V,L)
    return mu,ss,nsteps,probn

def optatBcfL(B,cfn,L,bnds=(-8,-2),ax=None):
    var = np.logspace(bnds[0],bnds[1],100)

    mu,ss =np.empty_like(var),np.empty_like(var)
    psteps = np.empty([3,var.size])

    for i,v in enumerate(var):
        c = cfn(v)
        mu[i],ss[i],nsteps,probn = perfatBCVL(B,c,v,L)
        psteps[:,i] = steps(nsteps,probn)

    iopt = np.argmin(np.exp(mu+0.5*ss))
    vopt = var[iopt]
    muopt = mu[iopt]
    ssopt = ss[iopt]
    stepsopt = psteps[:,iopt]
    if not ax is None:
        ax.plot(var,np.exp(mu),cols[0])
        ax.plot(var,np.exp(mu+2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(var,np.exp(mu-2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(var,np.exp(mu+0.5*ss),cols[0],linestyle='-.')
        ax.plot(vopt,np.exp(muopt+0.5*ssopt),color=cols[0],marker='o')
        axt = ax.twinx()
        axt.plot(var,psteps[1,:],cols[1])
        axt.plot(var,psteps[0,:],cols[1],linestyle='--',)
        axt.plot(var,psteps[2,:],cols[1],linestyle='--')
        ax.set_xscale('log')
        ax.set_yscale('log')
    return vopt,muopt,ssopt,stepsopt
#############################################################
fig,ax = plt.subplots(nrows=2,ncols=1)

B = 60*60
C = 60.
nsteps,probn = probsteps(B,C,ax=ax[0])
perfatl(nsteps,probn,1e-6,0.5,ax=ax[1])
fig.savefig('figs/perfplot.png')
plt.close(fig)
#############################################################
fig,ax = plt.subplots()
def cfn(v):
    return 60*1e-6/v

vopt,muopt,ssopt,stepsopt = optatBcfL(B,cfn,0.5,ax=ax)

fig.savefig('figs/rsprediction.png')
plt.close(fig)
#############################################################
fig,ax = plt.subplots(nrows=3,ncols=1,sharex=True)

L = np.logspace(np.log10(0.1),np.log10(1.5),150)
var = np.empty_like(L)
mu  = np.empty_like(L)
ss  = np.empty_like(L)
step= np.empty([3,L.size])

for i,ls in tqdm.tqdm(enumerate(L),total=L.size):
    var[i], mu[i], ss[i], step[:,i] = optatBcfL(B,cfn,ls)


ax[1].plot(L,var,cols[2])
ax[1].set_yscale('log')

ax[2].plot(L, step[1,:],cols[1])
ax[2].plot(L, step[0,:],cols[1],linestyle='--',)
ax[2].plot(L, step[2,:],cols[1],linestyle='--')

ax[0].plot(L,np.exp(mu),cols[0])
ax[0].plot(L,np.exp(mu+2*np.sqrt(ss)),cols[0],linestyle='--')
ax[0].plot(L,np.exp(mu-2*np.sqrt(ss)),cols[0],linestyle='--')
ax[0].plot(L,np.exp(mu+0.5*ss),cols[0],linestyle='-.')
ax[0].set_yscale('log')
ax[0].set_xscale('log')
fig.savefig('figs/lsprediction')