import numpy as np
import scipy as sp
from scipy import integrate
import copy
import gpbo
import time
import tqdm
import pandas as pd
from matplotlib import pyplot as plt
import pickle
ki = gpbo.core.GPdc.MAT52
Dim=2
lb = np.array([-1.]*Dim)
ub = np.array([1.]*Dim)

mpri = sp.array([2.]+[3.]*Dim)
spri = sp.array([0.5]+[0.15]*Dim)


def i3ntegrate23(f,bx,by,bz):
    L=[]
    def w(x,y,z):
        return [g(xi,y,z) for xi in x]
    def g(x,y,z):
        ll = f(x,y,z)
        L.append(ll)
        return np.exp(ll)
    def ix(y,z,o):
        I=[]
        for yi in y:
            I.append(integrate.fixed_quad(lambda x:w(x,yi,z),bx[0],bx[1],n=o)[0])
        return np.array(I)
    def ixy(z,o):
        I=[]
        for zi in z:
            I.append(integrate.fixed_quad(lambda y:ix(y,zi,o),by[0],by[1],n=o)[0])
        return np.array(I)
    def ixyz(o):
        I = integrate.fixed_quad(lambda z:ixy(z,o),bz[0],bz[1],n=o)
        return I[0]
    i2 = ixyz(2)
    i3 = ixyz(3)

    e=np.abs(i2-i3)
    lvm = max(L)
    return i3,e,lvm
class rect(object):
    def __init__(self,lb,ub,f,key=[]):
        self.key=key
        self.lb=lb
        self.ub=ub
        self.f=f
        i,e,l = i3ntegrate23(f,[lb[0],ub[0]],[lb[1],ub[1]],[lb[2],ub[2]])
        self.i=i
        self.e=e
        if e==0:
            self.e=l
        return
    def split(self):
        ax = np.argmax(np.abs(np.array(self.ub)-np.array(self.lb)))
        split = 0.5*(self.ub[ax]+self.lb[ax])
        lb0 = copy.copy(self.lb)
        lb1 = copy.copy(self.lb)
        lb1[ax]=split
        ub0 = copy.copy(self.ub)
        ub0[ax]=split
        ub1 = copy.copy(self.ub)
        return rect(lb0,ub0,self.f,key=self.key+[[ax,0]]),rect(lb1,ub1,self.f,key=self.key+[[ax,1]])
def opt(f):
    R=[rect([1e-3,1e-3,1e-3],[1e2,5.,5.],f)]
    A = copy.copy(R)
    for i in tqdm.tqdm(range(100)):
        s = np.argmax([r.e for r in R])
        rs = R.pop(s)
        r0,r1 = rs.split()
        R.append(r0)
        R.append(r1)
        A.append(r0)
        A.append(r1)
    I = np.sum([r.i for r in R])
    E = np.sum([max(0,r.e) for r in R])

    return I,E,[r.key for r in A]


#X = np.random.uniform(-1,1,[100,Dim])
#Y = sp.vstack([gpbo.core.objectives.shiftbraninojf(X[i,:],**{})[0] for i in range(X.shape[0])])
#S = np.ones_like(Y)*1e-6
#D = [[sp.NaN]]*Y.size

def readf(f):
    names = (open(f).readline().strip('\n')+''.join([',q{}'.format(i) for i in range(5)])).replace(' ','')
    df = pd.read_csv(f,names=names.split(','),skiprows=1,engine='c')
    n = len(df['y'].values)
    X = np.vstack([df['x{}'.format(i)].values[:n] for i in range(Dim)]).T
    Y = df['y'].values[:n].reshape([X.shape[0],1])
    S = df['s'].values[:n].reshape([X.shape[0],1])
    D = [[np.NaN]]*n
    return X,Y,S,D
X,Y,S,D = readf('results/pesfsb_0.csv')
global offset
offset=0
def getdivisions(n):
    #ql=lambda Q:-gpbo.core.GPdc.GP_LKonly(X[:n,:],Y[:n,:],S[:n,:],D[:n],gpbo.core.GPdc.kernel(ki,Dim,np.exp(np.array(Q)))).plk(mpri,spri,shape='gamma')
    #l_=np.inf
    #for i in range(50):
    #    hinit = np.log(np.array([np.random.uniform(0.1,100),np.random.uniform(0.05,1.5),np.random.uniform(0.05,1.5)]))
    #    l = ql(hinit)
     #   if l<l_:
     #       l_=l
     #       H=hinit
    #res=sp.optimize.minimize(ql,hinit,method='BFGS',options={'maxiter':100,'gtol':0.001})
    #print(res)
    #offset=res.fun
    H=[]
    P=[]
    def f(x,y,z):
        H.append([x,y,z])
        p= gpbo.core.GPdc.GP_LKonly(X[:n,:],Y[:n,:],S[:n,:],D[:n],gpbo.core.GPdc.kernel(ki,Dim,np.array([x,y,z]))).plk(mpri,spri,shape='gamma')
        #print(x,y,z,p,np.exp(p))
        P.append(p)
        global offset
        return p-offset
    I,E,K = opt(f)
    if np.isnan(I) or np.isinf(I):
        pass
    Ha = np.vstack(H)
    #print('mean {}\nmedian {}\nstd{}'.format(np.mean(Ha,axis=0),np.median(Ha,axis=0),np.sqrt(np.var(Ha,axis=0))))
    i = np.argmax(P)
    global offset
    print('offset {:.3g} maxp {:.3g} I {:.3g} E {:.3g}'.format(offset,P[i],I,E))

    offset = P[i]
    #print(H[i])
    #print(res.x)
    return K




def similarity(k0,K):
    S = np.zeros(len(k0))
    for kc in K:
        S = np.logical_or(S,np.array([k in kc for k in k0]))
    sim = np.sum(S)/float(S.size)
    if sim<0.5:
        pass
    return sim
#Kprev = getdivisions(10)
#for i in range(11,21):
#    K = getdivisions(i)
#    print(i,similarity(Kprev,K))
#    Kprev=K
def experiment():
    allK=[]
    allS=[]
    #header [n,unique,last1,last2,last3,new]
    Data = np.zeros([3,6])

    for e,i in enumerate(range(10,13)):
        Data[e,0]=i
        K0 = getdivisions(i)
        if len(allK)>=1:
            Data[e,2] = similarity(K0,allK[-1:])
        else:
            Data[e,2]=np.NaN
        if len(allK)>=2:
            Data[e,3] = similarity(K0,allK[-2:])
        else:
            Data[e,3]=np.NaN
        if len(allK)>=3:
            Data[e,4] = similarity(K0,allK[-3:])
        else:
            Data[e,4]=np.NaN
        Data[e,5] = similarity(K0,[allS])
        print('step {} prev1 {:.4g} prev2 {:.4g} prev3 {:.4g} ever {:.4g}'.format(i,similarity(K0,allK[-1:]),similarity(K0,allK[-2:]),similarity(K0,allK[-3:]),similarity(K0,[allS])))
        allK.append(K0)
        for k in K0:
            if not k in allS:
                allS.append(k)

        Data[e,1]=len(allS)
    return Data,K0

def plots(Data,K0):

    f,a = plt.subplots(nrows=1,ncols=1,sharex=True)
    a.plot(Data[:,0]+1,35*Data[:,1],'b',label='O(n^2) operations')
    a.plot(Data[:,0],35*len(K0)*(1-Data[:,5]),'r',label='O(n^3) operations')
    a.set_ylabel('Hyperparameter count')
    a.legend()
    f.savefig('figs/adaptive.png')

    f,a = plt.subplots(nrows=1,ncols=1,sharex=True)
    a.plot(Data[:,0],1-Data[:,5],'b',label='All previous steps')
    #a.plot(Data[:,0],1-Data[:,2],'r',label='1-step')
    a.plot(Data[:,0],1-Data[:,4],'r',label='3 previous steps')
    a.plot(Data[:,0],1-Data[:,2],'g',label='1 previous step')
    #a.set_yscale('symlog',linthreshy=1./len(K0))
    a.set_xlabel('Iteration')
    a.set_ylabel('Fraction repeated \nhyperparameter values')
    a.legend()
    f.savefig('figs/adaptivefrac.png')

if __name__=="__main__":
    if False:
        #run the experiment
        Data,K0 = experiment()
        pickle.dump([Data,K0],open('data/expdata.p','w'))

    plots(*pickle.load(open('data/expdata.p','r')))
