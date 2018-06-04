import numpy as np
import scipy as sp
from scipy.optimize import minimize_scalar
from matplotlib import pyplot as plt
import dill as pickle
import sys
import tqdm
import os
import time
import gpbo
from collections import defaultdict
cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

_,__,OP = pickle.load(open('/home/mark/gpbo/exps/predictive/overhead/overmodel.p'))
from gpbo.exps.predictive.overhead.overheads import cvimodel
def OMI(X):
    return cvimodel(X, OP)

_,PP = pickle.load(open('/home/mark/gpbo/exps/predictive/performance/model.p'))
from gpbo.exps.predictive.performance.model import perfmodel
def PM(X,S,L):
    Xa = np.ones([X.size,2])*np.log10(S)
    Xa[:,0] = X
    #Xa[:,1] = np.log10(S)
    #X = X.reshape(-1,1)
    #truemean,log10std = perfmodel(np.hstack([X,np.log10(S)*np.ones_like(X)]),L,PP)
    truemean,log10std = perfmodel(Xa,L,PP)
    natlogmu = np.log(truemean)
    natlogvar = (np.log(10)*log10std)**2
    return natlogmu,natlogvar



def probsteps(B,C,ax=None):
    def probstop(n):
        cm,cv = OMI(n)
        pneg = sp.stats.norm.cdf(0.,loc=cm,scale=np.sqrt(cv))
        pstop = (sp.stats.norm.cdf(np.maximum(0,B-C*(n+10)),loc=cm,scale=np.sqrt(cv))-pneg)/(1-pneg)
        return pstop

    nstepslog = np.logspace(1,np.log10(max(10,1+(B-10*C)/C)),250)
    probnlog = probstop(nstepslog)
    approxmaxarg = np.argmax(probnlog<0.125*1e-3)
    approxmax = nstepslog[approxmaxarg]
    nsteps = np.arange(int(approxmax))

    #nsteps = np.arange(max(1,min(200,1+(B-10*C)/C)))

    #bover = np.maximum(0,B-C*(nsteps+10))
    #cmean,cv = OMI(nsteps)
    #cstd = np.sqrt(cv)
    #cneg = sp.stats.norm.cdf(0.,loc=cmean,scale=cstd)
    cnstep = probstop(np.arange(int(approxmax)))
    last=1
    nstepthin=[]
    cnstepthin=[]
    for i in range(int(approxmax)):
        if cnstep[i]<last-0.25*1e-3:
            last = cnstep[i]
            nstepthin.append(i)
            cnstepthin.append(cnstep[i])
    cnstep = np.array(cnstepthin)
    nsteps = np.array(nstepthin)

    if len(nsteps)==0:
        return np.array([max(0,int(B/C))]),np.array([1])
    #(sp.stats.norm.cdf(bover,loc=cmean,scale=cstd)-cneg)/(1-cneg)

    cnstep[np.isnan(cnstep)]=1
    cnstep = np.minimum.accumulate(cnstep)

    pnstep = np.nan_to_num(np.hstack([-np.diff(cnstep),0]))
    if len(pnstep)==1:
        pnstep=np.ones(1)
    pnstep = pnstep/np.sum(pnstep)
    #initialization offset
    nsteps+=10
    if not ax is None:
        ns = np.linspace(0,nsteps[-1],150)
        cm,cv = OMI(ns)
        ax.plot(ns,cm,cols[0])
        ax.plot(ns,cm-2*np.sqrt(cv),cols[0],linestyle='--')
        ax.plot(ns,cm+2*np.sqrt(cv),cols[0],linestyle='--')

        ax.plot([0.,B/C],[B,0.],cols[0],linewidth=0.5)
        axt = ax.twinx()
        axt.plot(nsteps,cnstep,cols[1],linewidth=0.75)
    #if np.any(np.isnan(pnstep)):
    #    pass
    return nsteps, pnstep

def muss2mv(mu,ss):
    #natural log
    m = np.exp(mu+0.5*ss)
    v = (np.exp(ss)-1)*m**2#np.exp(2*mu+ss)
    return m,v
def mv2muss(m,v):
    #natural log
    #mu = np.log(m/np.sqrt(1+v/m**2))
    ss = np.log(1+v/m**2)
    mu = np.log(m)-0.5*ss
    return mu,ss

def marginalizelognorm(mu,ss,pn):
    #pn = p/np.sum(p)
    m,v = muss2mv(mu,ss)
    margm = np.sum(pn*m)
    margv = np.sum(pn*(m-margm)**2 + pn*v)
    margmu,margss = mv2muss(margm,margv)
    return margmu, margss

def perfatl(nsteps,probn,S,L,ax=None):
    #probn = np.zeros_like(probn)
    #probn[10]=1.
    xsteps = np.linspace(0,nsteps[-1],150)
    mu,ss = PM(xsteps,S,L)
    #lognormal mean variance in natural log

    mum,ssm = PM(nsteps,S,L)
    margmu,margss = marginalizelognorm(mum,ssm,probn)

    if not ax is None:
        ax.plot(xsteps,np.exp(mu),cols[0])
        ax.plot(xsteps,np.exp(mu+2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(xsteps,np.exp(mu-2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(xsteps,np.exp(mu+0.5*ss),cols[0],linestyle='-.')

        axt = ax.twinx()
        axt.plot(nsteps,probn/np.hstack([nsteps[0],np.diff(nsteps)]),cols[1],linewidth=0.75)

        ax.set_yscale('log')
        iplt = np.argmax(probn)
        xplt = nsteps[iplt]
        m,v = muss2mv(margmu,margss)
        ax.plot(xplt,m,'rx')
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
    #print(V,'steplen',len(nsteps))
    mu,ss = perfatl(nsteps,probn,V,L)
    return mu,ss,nsteps,probn

def perfatBCVoverL(B,C,V,L,pL):
    nsteps,probn = probsteps(B,C)
    #print(len(nsteps))
    n = pL.size
    mu = np.empty(n)
    ss = np.empty(n)
    for i in range(n):
        mu[i], ss[i] = perfatl(nsteps,probn,V,L[i])
    margmu, margss = marginalizelognorm(mu,ss,pL)
    return margmu, margss, nsteps, probn

def optatBcfL(B,cfn,L,bnds=(-6,-1),ax=None,axt=None):
    def minfn(v):
        c = cfn(10**v)
        mu,ss,nsteps,probn = perfatBCVL(B,c,10**v,L)
        m,v = muss2mv(mu,ss)
        return m
    print('optimize...')
    res = minimize_scalar(minfn, bounds=bnds, method='bounded',options={'maxiter':140})
    vopt = 10**res.x
    print('vopt = {}'.format(vopt))
    copt = cfn(vopt)
    muopt,ssopt,nsteps,probn = perfatBCVL(B,copt,vopt,L)
    m,v = muss2mv(muopt,ssopt)
    stepsopt = steps(nsteps,probn)
    Eover = B-np.sum(nsteps*probn*copt)

    if not ax is None:
        var = np.logspace(bnds[0],bnds[1],200)
        mu,ss =np.empty_like(var),np.empty_like(var)
        psteps = np.empty([3,var.size])
        for i,v in tqdm.tqdm(enumerate(var),total=len(var)):
            c = cfn(v)
            mu[i],ss[i],nsteps,probn = perfatBCVL(B,c,v,L)
            psteps[:,i] = steps(nsteps,probn)
        ax.plot(var,np.exp(mu),cols[0])
        ax.plot(var,np.exp(mu+2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(var,np.exp(mu-2*np.sqrt(ss)),cols[0],linestyle='--')
        ax.plot(var,np.exp(mu+0.5*ss),cols[0],linestyle='-.')
        ax.plot(vopt,np.exp(muopt+0.5*ssopt),color=cols[0],marker='o',linestyle='None')
        #axt = ax.twinx()
        axt.plot(var,psteps[1,:],cols[1])
        axt.plot(var,psteps[0,:],cols[1],linestyle='--',)
        axt.plot(var,psteps[2,:],cols[1],linestyle='--')
        axt.grid('False')
        ax.set_xscale('log')
        ax.set_yscale('log')

    R = {'mu':muopt, 'ss':ssopt, 'Eover':Eover, 'Esteps':stepsopt[1], 'B':B, 'c':copt, 'obsvar':vopt,'Rmean':m,'Rvar':v}
    return R

def optatBcfoverL(B,cfn,L,pL,bnds=(-6,-1),ax=None,axt=None):
    sys.stdout.flush()
    #with tqdm.tqdm() as pbar:
    def minfn(v):
 #       pbar.update()
        c = cfn(10**v)
        mu,ss,nsteps,probn = perfatBCVoverL(B,c,10**v,L,pL)
        m,v = muss2mv(mu,ss)
        return m
    res = minimize_scalar(minfn, bounds=bnds, method='bounded',options={'maxiter':40})
    vopt = 10**res.x
    sys.stdout.flush()
    #print('vopt = {}\n'.format(vopt))
    copt = cfn(vopt)
    muopt,ssopt,nsteps,probn = perfatBCVoverL(B,copt,vopt,L,pL)
    m,v = muss2mv(muopt,ssopt)
    stepsopt = steps(nsteps,probn)
    Eover = B-np.sum(nsteps*probn*copt)
    if not ax is None:
        var = np.logspace(bnds[0],bnds[1],200)
        mu,ss =np.empty_like(var),np.empty_like(var)
        psteps = np.empty([3,var.size])
        for i,v in tqdm.tqdm(enumerate(var),total=len(var)):
            c = cfn(v)
            mu[i],ss[i],nsteps,probn = perfatBCVoverL(B,c,v,L,pL)
            #ovm,ovv = OMI(nsteps-10)
            #meanover = np.sum(probn*ovm)
            #varover = np.sum(probn*ovv)+np.sum(probn*(ovm-meanover)**2)
            #psteps[:,i] = steps(nsteps,probn)
            alow = np.argmax(np.cumsum(probn)>0.025)
            ahigh = np.argmax(np.cumsum(probn)>0.975)
            evc = min(B,np.sum(c*nsteps*probn))
            psteps[1,i] = (B-evc)/float(B)
            evl = min(B,np.sum(c*nsteps[alow]))
            evh = min(B,np.sum(c*nsteps[ahigh]))
            psteps[0,i] = (B-evl)/float(B)
            psteps[2,i] = (B-evh)/float(B)
        ax.plot(var,np.exp(mu),cols[0],label=' Mean $\\log_{10}$ Regret Prediction')
        ax.fill_between(var,np.exp(mu-2*np.sqrt(ss)),np.exp(mu+2*np.sqrt(ss)),facecolor=cols[0],edgecolor='None',alpha=0.2)
        ax.plot(var,np.exp(mu+0.5*ss),cols[0],linestyle='-.',label='Mean Regret Prediction')
        ax.plot(vopt,np.exp(muopt+0.5*ssopt),color=cols[0],marker='o',linestyle=None)
        #axt = ax.twinx()
        am = np.argmax(vopt<var)
        axt.plot(var,psteps[1,:],cols[1],label='Overhead Fraction Prediction')
        axt.fill_between(var,psteps[0,:],psteps[2,:],facecolor=cols[1],edgecolor='None',alpha=0.2)
        axt.plot(var[am],psteps[1,am],color=cols[1],marker='o')
        axt.set_ylim(0,1)
        #print(var,np.exp(mu))

    R = {'mu':muopt, 'ss':ssopt, 'Eover':Eover, 'Esteps':stepsopt[1], 'B':B, 'c':copt, 'obsvar':vopt,'Rmean':m,'Rvar':v}
    return R

def plottruedata(path,base,axR,axN):

    names = [f for f in os.listdir(path) if f.startswith(base)]
    D = [gpbo.optimize.readoptdata(os.path.join(path,f)) for f in names]
    for optdata in D:
        N = len(optdata['index'])
        R = optdata['trueyatxrecc'].values[-1]
        n = optdata['s'][0]
        axR.plot(n,R,'b.')
        axN.plot(n,N,'r.')
    return

def plotmargdata(path,base,ls,lw,axR,axN):
    result = defaultdict(lambda :[])
    for i in tqdm.tqdm(range(len(ls))):
        print(i)
        names = [f for f in os.listdir(path) if f.startswith(base+str(ls[i]))]
        D = [gpbo.optimize.readoptdata(os.path.join(path,f)) for f in names]
        for optdata in D:
            N = len(optdata['index'])
            R = optdata['trueyatxrecc'].values[-1]
            n = optdata['s'][0]
            result[n].append([R,N,lw[i]])
    for k in result.keys():
        data = np.vstack(result[k])
        mean = np.sum(data[:,0]*data[:,2]/np.sum(data[:,2]))
        numstep = np.sum(data[:,1]*data[:,2]/np.sum(data[:,2]))
        axR.plot(k,mean,'b.')
        axN.plot(k,numstep,'r.')
    return

if __name__ == "__main__":
    #############################################################
    fig,ax = plt.subplots(nrows=2,ncols=1,sharex=True)
    B = 10*3600
    C = 60.
    nsteps,probn = probsteps(B,C,ax=ax[0])
    #print(len(nsteps))
    #print(probsteps2(60*60,60)[0])
    #print(probsteps2(60*60,6)[0])
    #print(probsteps2(3*60*60,6)[0])


    #sys.exit(0)
    perfatl(nsteps,probn,1e-4,0.5,ax=ax[1])
    fig.savefig('figs/perfplot.png')
    plt.close(fig)
    #############################################################
    print('rsplot')
    fig,ax = plt.subplots(nrows=4,ncols=1,figsize=[5,8],sharex=True)
    def cfn(v):
        #return 60*1e-6/v
        return B*0.01*1e-4/v
    axt = [a.twinx() for a in ax]
    R = optatBcfL(B,cfn,0.1,ax=ax[0],axt=axt[0],bnds=(-8,-2))
    R = optatBcfL(B,cfn,0.5,ax=ax[1],axt=axt[1],bnds=(-8,-2))
    R = optatBcfL(B,cfn,1.2,ax=ax[2],axt=axt[2],bnds=(-8,-2))

    #plottruedata('results','eihyp_3_100_',ax[0],axt[0])
    #plottruedata('results','eihyp_3_500_',ax[1],axt[1])
    #plottruedata('results','eihyp_3_1200_',ax[2],axt[2])
    #lset = np.array([100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400])
    #wts = sp.stats.gamma.pdf(lset/1000.,3.,scale =0.2)
    #wts = wts/np.sum(wts)

    #plotmargdata('margresults','eihyp_3_',lset,wts,ax[3],axt[3])
    #L = sp.stats.gamma.rvs(3,scale=0.15,size=500)
    L = sp.stats.gamma.ppf(np.linspace(0,1,1002)[1:-1],4.,scale =0.2)
    p = np.ones_like(L)/float(L.size)

    R = optatBcfoverL(B,cfn,L,p,ax=ax[3],axt=axt[3],bnds=(-8,-2))

    #print(R)
    fig.savefig('figs/rsprediction.png')
    plt.close(fig)


    sys.exit(0)
    #############################################################
    fig,ax = plt.subplots(nrows=3,ncols=1,sharex=True)

    L = np.logspace(np.log10(0.1),np.log10(1.5),150)
    var = np.empty_like(L)
    mu  = np.empty_like(L)
    ss  = np.empty_like(L)
    step= np.empty([3,L.size])

    for i,ls in tqdm.tqdm(enumerate(L),total=L.size):
        var[i], mu[i], ss[i], step[:,i] = optatBcfL(B,cfn,ls)

    #overfrac = np.array([np.sum(OM(np.arange(s))[0]) for s in step[1,:]])/float(B)


    ax[1].plot(L,var,cols[2])
    ax[1].set_yscale('log')

    #ax[1].twinx().plot(L,overfrac,cols[3])
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