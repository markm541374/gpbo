import numpy as np
import scipy as sp
import gpbo
import copy
import tqdm
import pandas as pd
import cho_update as cu
import time
from matplotlib import pyplot as plt
import pickle
def logintegrator2(f,lb,ub):
    #2 order integral that uses log values
    r=ub[0]-lb[0]
    lr = np.log(r)
    if len(lb)==1:
        il = f(lb[0])
        iu = f(ub[0])
    else:
        il = logintegrator2(lambda *args:f(lb[0],*args),lb[1:],ub[1:])
        iu = logintegrator2(lambda *args:f(ub[0],*args),lb[1:],ub[1:])
    mx = max(il,iu)
    return np.log(0.5*np.exp(il-mx)+0.5*np.exp(iu-mx))+mx+lr

def logintegrator3(f,lb,ub):
    #2 order integral that uses log values
    r=ub[0]-lb[0]
    lr = np.log(r)
    if len(lb)==1:
        il = f(lb[0])
        im = f(0.5*(lb[0]+ub[0]))
        iu = f(ub[0])
    else:
        il = logintegrator3(lambda *args:f(lb[0],*args),lb[1:],ub[1:])
        im = logintegrator3(lambda *args:f(0.5*(lb[0]+ub[0]),*args),lb[1:],ub[1:])
        iu = logintegrator3(lambda *args:f(ub[0],*args),lb[1:],ub[1:])
    mx = max(il,im,iu)
    return np.log(0.25*np.exp(il-mx)+0.5*np.exp(im-mx)+0.25*np.exp(iu-mx))+mx+lr

class rect(object):
    def __init__(self,lb,ub,intfn,key=[]):
        self.key=key
        self.lb=lb
        self.ub=ub
        self.intfn=intfn
        self.depth=len(key)
        lf2 = logintegrator2(intfn,lb,ub)
        lf3 = logintegrator3(intfn,lb,ub)
        self.i=lf3
        mx = max(lf2,lf3)
        self.e = np.log(np.abs(np.exp(lf3-mx)-np.exp(lf2-mx)))+mx
        return
    def split(self):
        ax = self.depth%len(self.ub)#np.argmax(np.abs(np.array(self.ub)-np.array(self.lb)))
        split = 0.5*(self.ub[ax]+self.lb[ax])
        lb0 = copy.copy(self.lb)
        lb1 = copy.copy(self.lb)
        lb1[ax]=split
        ub0 = copy.copy(self.ub)
        ub0[ax]=split
        ub1 = copy.copy(self.ub)
        return rect(lb0,ub0,self.intfn,key=self.key+[[ax,0]]),rect(lb1,ub1,self.intfn,key=self.key+[[ax,1]])

def opt(f,lb,ub):
    R=[rect(lb,ub,f)]
    A = copy.copy(R)
    for i in tqdm.tqdm(range(200)):
        s = np.argmax([r.e for r in R])
        rs = R.pop(s)
        r0,r1 = rs.split()
        R.append(r0)
        R.append(r1)
        A.append(r0)
        A.append(r1)
    LI = np.array([r.i for r in R])
    mx = LI.max()
    I = np.sum(np.exp(LI-mx))*np.exp(mx)

    EI = np.array([r.e for r in R])
    mx = EI.max()
    E = np.sum(np.exp(EI-mx))*np.exp(mx)

    return I,E,[r.key for r in A]

def runadaptive(X,Y,n):
    Dim=X.shape[1]
    LL = cu.gp_llk(copy.deepcopy(X),copy.deepcopy(Y),gpbo.core.GPdc.MAT52,Dim,1e-6)
    prev=0

    totalpts=[]
    newpts=[]
    effectivepts=[]
    activepts=[]
    npts=[]
    meanupdates=[]
    updatetime=[]
    touchedpts=[]
    for i in range(10,n):
        npts.append(i)
        print('i={}'.format(i))
        H=[]
        L={}
        U=[]
        def f(a,l1,l2):
            H.append([a,l1,l2])
            l,u= LL[(a,l1,l2)].llk(i)
            L[(a,l1,l2)]=l
            if u>0:
                U.append(u)
            return l
        t0=time.clock()
        Q = opt(f,[1e-6,1e-6,1e-6],[100.,4.,4.])
        t1=time.clock()
        updatetime.append(t1-t0)
        hyp = np.vstack(H)
        print('mean {}'.format(hyp.mean(axis=0)))
        print('std  {}'.format(np.sqrt(hyp.var(axis=0))))

        totalpts.append(len(LL))
        newpts.append(len(LL)-prev)
        effectivepts.append(np.sum(np.exp(np.array(L.values()))>0.01*np.exp(max(L.values()))))
        activepts.append(len(L))
        meanupdates.append(np.mean(U))
        touchedpts.append(len(U))
        print(totalpts[-1],newpts[-1],newpts[-1]/float(totalpts[-1]))
        print(activepts[-1],effectivepts[-1],meanupdates[-1])
        prev=totalpts[-1]
    d = {'npts':npts,
         'newpts':newpts,
         'effectivepts':effectivepts,
         'activepts':activepts,
         'touchedpts':touchedpts,
         'meanupdates':meanupdates,
         'updatetime':updatetime,
         'totalpts':totalpts}
    return d

def runfixg(X,Y,n,N):
    Dim=X.shape[1]
    LL = cu.gp_llk(copy.deepcopy(X),copy.deepcopy(Y),gpbo.core.GPdc.MAT52,Dim,1e-6)
    prev=0
    gi = [np.logspace(-2,2,N)]*(Dim+1)
    G = np.vstack(np.meshgrid(*gi)).reshape(Dim+1,-1).T
    totalpts=[]
    newpts=[]
    effectivepts=[]
    activepts=[]
    npts=[]
    meanupdates=[]
    updatetime=[]
    touchedpts=[]
    for i in range(10,n):
        npts.append(i)
        print('i={}'.format(i))
        H=[]
        L={}
        U=[]
        def f(a,l1,l2):
            H.append([a,l1,l2])
            l,u= LL[(a,l1,l2)].llk(i)
            L[(a,l1,l2)]=l
            if u>0:
                U.append(u)
            return l
        t0=time.clock()
        for j in range(G.shape[0]):
            loglike=f(*[k for k in G[j,:]])
        t1=time.clock()
        updatetime.append(t1-t0)
        hyp = np.vstack(H)
        print('mean {}'.format(hyp.mean(axis=0)))
        print('std  {}'.format(np.sqrt(hyp.var(axis=0))))

        totalpts.append(len(LL))
        newpts.append(len(LL)-prev)
        effectivepts.append(np.sum(np.exp(np.array(L.values()))>0.01*np.exp(max(L.values()))))
        activepts.append(len(L))
        meanupdates.append(np.mean(U))
        touchedpts.append(len(U))
        print(totalpts[-1],newpts[-1],newpts[-1]/float(totalpts[-1]))
        print(activepts[-1],effectivepts[-1],meanupdates[-1])
        prev=totalpts[-1]

    d = {'npts':npts,
         'newpts':newpts,
         'effectivepts':effectivepts,
         'touchedpts':touchedpts,
         'activepts':activepts,
         'meanupdates':meanupdates,
         'updatetime':updatetime,
         'totalpts':totalpts}
    return d

def runslice(X,Y,n,N):
    Dim=X.shape[1]
    prev=0
    totalpts=[]
    newpts=[]
    effectivepts=[]
    activepts=[]
    npts=[]
    meanupdates=[]
    updatetime=[]
    touchedpts=[]
    for i in range(10,n):
        #LL = cu.gp_llk(copy.deepcopy(X),copy.deepcopy(Y),gpbo.core.GPdc.MAT52,Dim,1e-6)
        npts.append(i)
        print('i={}'.format(i))
        H=[]
        L={}
        U=[]
        def f(hyp):
            a,l1,l2=[np.exp(h) for h in hyp]
            if np.any(np.exp(hyp)>[100,4,4]) or np.any(np.exp(hyp)<[1e-6,1e-6,1e-6]):
                return -1e99
            H.append([a,l1,l2])
            k=gpbo.core.GPdc.kernel(gpbo.core.GPdc.MAT52,Dim,np.array((a,l1,l2)))
            l,u=cu.lazyllk(X,Y,k,i,1e-6).llk(i)
            if u>0:
                U.append(u)
            return l
        t0=time.clock()
        start=np.array([1.,0.,0.])
        W = gpbo.core.slice.slice_sample(f,start,N,[0.1]*(Dim+1),burn=20,subsam=1)
        t1=time.clock()
        LK = [f(W[j,:]) for j in range(W.shape[0])]
        updatetime.append(t1-t0)
        try:
            hyp = np.vstack(H)
            print('mean {}'.format(hyp.mean(axis=0)))
            print('std  {}'.format(np.sqrt(hyp.var(axis=0))))
        except:
            print('hyp invalid {} {}'.format(W,H))
        totalpts.append(len(H))
        effectivepts.append(np.sum(np.exp(np.array(LK))>0.01*np.exp(max(LK))))
        activepts.append(W.shape[0])
        meanupdates.append(i)
        touchedpts.append(len(U))
        newpts.append(len(U))
        print(totalpts[-1],npts[-1],1)
        print(W.shape[0],effectivepts[-1],meanupdates[-1])
    d = {'npts':npts,
         'newpts':newpts,
         'effectivepts':effectivepts,
         'activepts':activepts,
         'touchedpts':touchedpts,
         'meanupdates':meanupdates,
         'updatetime':updatetime,
         'totalpts':totalpts}
    return d
def readf(f,Dim):
    names = (open(f).readline().strip('\n')+''.join([',q{}'.format(i) for i in range(5)])).replace(' ','')
    df = pd.read_csv(f,names=names.split(','),skiprows=1,engine='c')
    n = len(df['y'].values)
    X = np.vstack([df['x{}'.format(i)].values[:n] for i in range(Dim)]).T
    X = np.ascontiguousarray(X)
    Y = df['y'].values[:n].reshape([X.shape[0],1])
    S = df['s'].values[:n].reshape([X.shape[0],1])
    D = [[np.NaN]]*n
    return X,Y,S,D
X,Y,_,__ = readf('results/pesfsb_0.csv',2)

if False:
    d0 = runadaptive(X,Y,100)
    pickle.dump([d0],open('data/adapt.p','w'))
    d1 = runfixg(X,Y,100,14)
    pickle.dump([d1],open('data/grid.p','w'))
    d2 = runslice(X,Y,100,200)
    pickle.dump([d2],open('data/slice.p','w'))
else:
    d0 = pickle.load(open('dataheron/adapt.p','r'))[0]
    d1 = pickle.load(open('dataheron/grid.p','r'))[0]
    d2 = pickle.load(open('dataheron/slice.p','r'))[0]

cmap = plt.rcParams['axes.prop_cycle'].by_key()['color']
d0['name']='Adaptive Quad'
d1['name']='Regular Grid'
d2['name']='Slice Sample'
d0['c']=cmap[0]
d1['c']=cmap[1]
d2['c']=cmap[2]

def plottime(D):
    print('plotting times')
    f,a = plt.subplots(1)
    [a.plot(d['npts'],d['updatetime'],color=d['c'],label=d['name']) for d in D]
    a.set_ylabel('Marginalization Time (s)')
    a.set_xlabel('Iteration')
    a.legend()
    #a.set_yscale('log')
    f.savefig('figs/overheadtime.pdf')

def plotops(D):
    print('plotting ops')
    f,a = plt.subplots(1)
    [a.plot(d['npts'],np.array(d['touchedpts'])-np.array(d['newpts']),color=d['c'],label=d['name']+' $\mathcal{O}(n^2)$ Updates',linestyle='--') for d in D]
    [a.plot(d['npts'],np.array(d['newpts']),color=d['c'],label=d['name']+' $\mathcal{O}(n^3)$ Evaluations') for d in D]
    a.set_ylabel('Hyperparameter LLK Evaluations')
    a.set_xlabel('Iteration')
    a.set_ylim(0,2*1e4)
    a.legend()
    #a.set_yscale('symlog',linthreshy=10000)
    f.savefig('figs/opcount.pdf')
def plotuseful(D):
    print('plotting useful')
    f,a = plt.subplots(1)
    [a.plot(d['npts'],d['activepts'],color=d['c'],label=d['name']+' Active') for d in D]
    [a.plot(d['npts'],d['effectivepts'],color=d['c'],label=d['name']+' Effective',linestyle='--') for d in D]
    a.set_ylabel('Hyperparameters Used (s)')
    a.set_xlabel('Iteration')
    a.legend(loc=1)
    f.savefig('figs/useful.pdf')

def plotmeanupdate(d):
    print('plotting meanupdate')
    f,a = plt.subplots(1)
    a.plot(d['npts'],d['meanupdates'],color=d['c'])
    a.set_xlabel('Iteration')
    a.set_ylabel('Mean Update Order')
    f.savefig('figs/updateorder.pdf')

#plottime([d0,d1,d2])
plotops([d0,d1,d2])
#plotuseful([d0,d1,d2])
#plotmeanupdate(d0)