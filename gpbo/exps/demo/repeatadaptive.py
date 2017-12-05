import numpy as np
import scipy as sp
from scipy import integrate
import copy
import gpbo
import time
import tqdm
from matplotlib import pyplot as plt
ki = gpbo.core.GPdc.MAT52
Dim=2
lb = np.array([-1.]*Dim)
ub = np.array([1.]*Dim)

mpri = sp.array([2.]+[3.]*Dim)
spri = sp.array([0.5]+[0.15]*Dim)


def i3ntegrate23(f,bx,by,bz):
    w = lambda x,y,z:[f(xi,y,z) for xi in x]
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
    #print(i3,i4)
    e=np.abs(i2-i3)
    return i3,e
class rect(object):
    def __init__(self,lb,ub,f,key=[]):
        self.key=key
        self.lb=lb
        self.ub=ub
        self.f=f
        i,e = i3ntegrate23(f,[lb[0],ub[0]],[lb[1],ub[1]],[lb[2],ub[2]])
        self.i=i
        self.e=e
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
    for i in tqdm.tqdm(range(1000)):
        s = np.argmax([r.e for r in R])
        rs = R.pop(s)
        r0,r1 = rs.split()
        R.append(r0)
        R.append(r1)
        A.append(r0)
        A.append(r1)
    I = np.sum([r.i for r in R])
    E = np.sum([r.e for r in R])
    #print([r.i for r in R])
    return I,E

def timing(n):
    X = np.random.uniform(-1,1,[n,Dim])
    Y = sp.vstack([gpbo.core.objectives.shiftbraninojf(X[i,:],**{})[0] for i in range(X.shape[0])])
    S = np.ones_like(Y)
    D = [[sp.NaN]]*Y.size
    H=[]
    P=[]
    def f(x,y,z):
        H.append([x,y,z])
        p= gpbo.core.GPdc.GP_LKonly(X,Y,S,D,gpbo.core.GPdc.kernel(ki,Dim,np.array([x,y,z]))).plk(mpri,spri,shape='gamma')
        #print(x,y,z,p,np.exp(p))
        P.append(p)
        return np.exp(p)
    print(opt(f))
    Ha = np.vstack(H)
    print('mean {}\nmedian {}\nstd{}'.format(np.mean(Ha,axis=0),np.median(Ha,axis=0),np.sqrt(np.var(Ha,axis=0))))
    i = np.argmax(P)
    print(H[i])

    #G = gpbo.core.PES.makeG(X,Y,S,D,ki, mpri,spri,50,prior='gamma')
    res=sp.optimize.minimize(lambda Q:-gpbo.core.GPdc.GP_LKonly(X,Y,S,D,gpbo.core.GPdc.kernel(ki,Dim,np.array(Q))).plk(mpri,spri,shape='gamma'),[10.,0.5,0.5],method='Nelder-Mead')
    print(res.x)
timing(80)


