import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from gpbo.core import objectives
import gpbo
import time
import pandas as pd
import sys

fout = 'results/b2support.csv'
fin = 'results/pesfsb_0.csv'

ki = gpbo.core.GPdc.MAT52
Dim=2
lb = np.array([-1.]*Dim)
ub = np.array([1.]*Dim)

mpri = np.array([1.]+[0.]*Dim)
spri = np.array([1.]*(Dim+1))
#mpri=np.array([2.]+[3.]*Dim)
#spri=np.array([0.5]+[0.15]*Dim)
def measures(n,f):
    print('measures with n={}'.format(n))
    names = (open(f).readline().strip('\n')+''.join([',q{}'.format(i) for i in range(5)])).replace(' ','')
    df = pd.read_csv(f,names=names.split(','),skiprows=1,engine='c')
    X = np.vstack([df['x{}'.format(i)].values[:n] for i in range(Dim)]).T
    Y = df['y'].values[:n].reshape([X.shape[0],1])
    S = df['s'].values[:n].reshape([X.shape[0],1])
    D = [[np.NaN]]*n

    #MAP = gpbo.core.GPdc.searchMAPhyp(np.copy(X), np.copy(Y), np.copy(S), D, mpri, spri, ki)
    #H = gpbo.core.ESutils.drawhyp_plk(X,Y,S,D,ki,mpri,spri,28,chains=1,prior='gamma')
    G = gpbo.core.PES.makeG(X,Y,S,D,ki, mpri,spri,28,prior='lognorm')

    hyp = sp.array([k.hyp for k in G.kf])
    hmean = sp.mean(hyp, axis=0)
    hstd = sp.sqrt(sp.var(hyp, axis=0))
    hmin = hyp.min(axis=0)
    hmax = hyp.max(axis=0)
    hmed = sp.median(hyp,axis=0)
    print('hyperparameters:\nmean {}\nmedian {}\nstd {}\nmin {}\nmax {}'.format(hmean,hmed,hstd,hmin,hmax))


    def kl(Z,G):
        m = Z.shape[0]
        c = np.zeros(m)
        M = np.zeros(10)
        for i in range(10):
            P = G.draw_post(Z,[[np.NaN]]*m,1000)
            A = np.argmin(P,axis=1)
            M[i]=np.mean(P.min(axis=1))
            for a in A:
                c[a]+=1
        p = (c+1)/(m+np.sum(c))
        q = 1./float(m)
        return np.sum(-p*np.log(q/p)),np.mean(M),np.sum(c==0)

    plot=False
    if plot:
        f,a = plt.subplots(nrows=4,ncols=2,figsize=[12,17])
    m = 1000
    def plotdraw(Z,a):
        a.plot(Z[:,0],Z[:,1],'b.')
        P = G.draw_post(Z,[[np.NaN]]*m,1000)
        A = np.argmin(P,axis=1)
        a.plot(Z[:,0][A],Z[:,1][A],'r.')
        return
    try:
        t0=time.clock()
        Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_LAPAPROT,para=12,pad_unif=False)
        T = time.clock()-t0
        laprot = list(kl(Z,G))
        laprot.append(T)
        if plot:
            plotdraw(Z,a[0,0])
    except:
        laprot = [np.Inf,np.Inf,np.Inf,np.Inf]
    try:
        t0=time.clock()
        Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_LAPAPROT,para=12,pad_unif=False,weighted=2)
        T = time.clock()-t0
        laprotw = list(kl(Z,G))
        laprotw.append(T)
        if plot:
            plotdraw(Z,a[0,1])
    except:
        laprotw = [np.Inf,np.Inf,np.Inf,np.Inf]

    try:
        t0=time.clock()
        Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_SLICEEI)
        T = time.clock()-t0
        ei = list(kl(Z,G))
        ei.append(T)
        if plot:
            plotdraw(Z,a[1,0])
    except:
        ei = [np.Inf,np.Inf,np.Inf,np.Inf]
    try:
        t0=time.clock()
        Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_SLICELCB,para=2)
        T = time.clock()-t0
        lcb = list(kl(Z,G))
        lcb.append(T)
        if plot:
            plotdraw(Z,a[1,1])
    except:
        lcb = [np.Inf,np.Inf,np.Inf,np.Inf]

    try:
        t0=time.clock()
        Z = gpbo.core.ESutils.draw_support(G.D, lb, ub, m,gpbo.core.ESutils.SUPPORT_UNIFORM)
        T = time.clock()-t0
        unif = list(kl(Z,G))
        unif.append(T)
        if plot:
            plotdraw(Z,a[2,0])
    except:
        unif = [np.Inf,np.Inf,np.Inf,np.Inf]

    if plot:
        l = 200
        x_ = sp.linspace(-1,1,l)
        y_ = sp.linspace(-1,1,l)
        z_ = sp.empty([l,l])
        s_ = sp.empty([l,l])
        for i in range(l):
            for j in range(l):
                m_,v_ = G.infer_diag_post(sp.array([y_[j],x_[i]]),[[sp.NaN]])
                z_[i,j] = m_[0,0]
                s_[i,j] = sp.sqrt(v_[0,0])
        CS = a[2,1].contour(x_,y_,z_,20)
        a[2,1].clabel(CS, inline=1, fontsize=10)
        CS = a[3,1].contour(x_,y_,s_,20)
        a[3,1].clabel(CS, inline=1, fontsize=10)
        f.savefig('dbout/h{}.svg'.format(n))
    plt.close(f)

    return laprot,laprotw,ei,lcb,unif


with open(fout,'w') as f:
    f.write('n,LRkl,LRy,LRm,LRt,LHkl,LHy,LHm,LHt,EIkl,EIy,EIm,EIt,CBkl,CBy,CBm,CBt,Ukl,Uy,Um,Ut\n')
for i in np.arange(10,80):
    with open(fout,'a') as f:
        res=np.hstack([[i],np.hstack(measures(i,fin))])
        f.write(','.join([str(j) for j in res])+'\n')

sys.exit(0)

npts = np.linspace(10,120,30).astype(int)
nv = npts.size
f,a = plt.subplots(nrows=3,ncols=1)
a[0].set_yscale('log')
for i in range(nv):
    laprot,laprotw,ei,lcb,unif = measures(npts[i])
    print(laprot,laprotw,ei,lcb,unif)

    a[0].plot(npts[i],laprot[0],'b.')
    a[0].plot(npts[i],laprotw[0],'r.')
    a[0].plot(npts[i],ei[0],'g.')
    a[0].plot(npts[i],lcb[0],'k.')
    a[0].plot(npts[i],unif[0],'m.')

    a[1].plot(npts[i],laprot[1],'b.')
    a[1].plot(npts[i],laprotw[1],'r.')
    a[1].plot(npts[i],ei[1],'g.')
    a[1].plot(npts[i],lcb[1],'k.')
    a[1].plot(npts[i],unif[1],'m.')

    a[2].plot(npts[i],laprot[2],'b.')
    a[2].plot(npts[i],laprotw[2],'r.')
    a[2].plot(npts[i],ei[2],'g.')
    a[2].plot(npts[i],lcb[2],'k.')
    a[2].plot(npts[i],unif[2],'m.')
    f.savefig('figs/qualityh.svg')