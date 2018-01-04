import numpy as np
import scipy as sp
from scipy.linalg import cho_factor
from scipy.linalg import cho_solve
import time
import tqdm
from matplotlib import pyplot as plt
from collections import defaultdict
import gpbo


class kupdate(object):
    def __init__(self,maxsize,x,y,k,jit=1e-6):
        self.k = k
        self.jit=jit
        self.C = np.zeros([maxsize,maxsize])
        self.X = np.zeros([maxsize,x.shape[1]])
        self.Y = np.zeros(maxsize)

        #self.M = np.zeros([maxsize,maxsize])
        self.count=0
        for i in range(x.shape[0]):
            self.update(x[i,:],y[i])
        return

    def update(self,x,y):
        self.X[self.count,:]=x
        self.Y[self.count]=y
        for i in range(self.count):
            #self.M[self.count,i]=self.k(x,self.X[i,:])
            self.C[self.count,i]=(self.k(x,self.X[i,:])-self.C[self.count,:i].dot(self.C[i,:i]))/self.C[i,i]

        #self.M[self.count,self.count]=self.k(x,x)+self.jit
        self.C[self.count,self.count] = self.k(x,x)+self.jit-np.sum(self.C[self.count,:self.count]**2)
        if self.C[self.count,self.count]<0.:
            raise ValueError
        self.C[self.count,self.count]=np.sqrt(self.C[self.count,self.count])

        self.count+=1

        #self.H = cho_factor(self.M[:self.count,:self.count],lower=True)
        #print(np.allclose(self.C[:self.count,:self.count],np.array(self.H[0])))
        #self.llk()
        return
    def llk(self,n):
        ym=np.mean(self.Y[:n])

        t2 = -0.5*(self.Y[:n]-ym).dot(cho_solve((self.C[:n,:n],True),self.Y[:n]-ym))
        t0 = -0.5*n*np.log(2*sp.pi)
        t1 = -np.sum(np.log(np.diag(self.C[:n,:n])))
        lk=t0+t1+t2
#        p= gpbo.core.GPdc.GP_LKonly(self.X[:n,:],self.Y[:n],self.jit*np.ones([n,1]),[[np.NaN]]*self.count,k).llk()
 #       print(p,lk)
        return t0+t1+t2

class lazyllk(object):
    def __init__(self,X,Y,k,maxsize,jit):
        self.X=X
        self.Y=Y
        self.cu=kupdate(maxsize,X[:1,:],Y[:1,:],k,jit=jit)
        self.cache={}
        return

    def llk(self,n):
        if n in self.cache.keys():
            return self.cache[n],0
        nupdates=n-self.cu.count
        while self.cu.count<n:
            self.cu.update(self.X[self.cu.count,:],self.Y[self.cu.count,:])
        L = self.cu.llk(n)
        self.cache[n]=L
        return L,nupdates


class keydefaultdict(defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret

def gp_llk(X,Y,ki,Dim,jit):

    def new(h):
        #print(h)
        k=gpbo.core.GPdc.kernel(ki,Dim,np.array(h))
        return lazyllk(X,Y,k,200,jit)

    return keydefaultdict(lambda h:new(h))


def main():
    ki = gpbo.core.GPdc.MAT52
    Dim=2
    k=gpbo.core.GPdc.kernel(ki,Dim,np.array([1.1,0.2,0.3]))
    X = np.random.uniform(-1,1,[400,Dim])
    Y = sp.vstack([gpbo.core.objectives.shiftbraninojf(X[i,:],**{})[0] for i in range(X.shape[0])])
    S = np.ones_like(Y)*1e-6
    D = [[sp.NaN]]*Y.size

    #def new(h):
    #    k=gpbo.core.GPdc.kernel(ki,Dim,np.array(h))
     #   return lazyllk(X,Y,k,200,1e-6)
    #G = kupdate(400,X[:0,:],Y[:0,0],k,jit=S[0,0])
    T=[]
    L = gp_llk(X,Y,ki,Dim,1e-6)
        #keydefaultdict(lambda h:new(h))
    _ = L[(1.1,0.2,0.3)]
    for i in tqdm.tqdm(range(100)):
        #if i==50:
        #    _=L[(1.1,0.2,0.3)].llk(70)
        t0=time.clock()
        #G.update(X[i,:],Y[i,0])
        l1=L[(1.1,0.2,0.3)].llk(i+1)[0]
        t1=time.clock()
    #    Gf = lakupdate(i+1,X[:(i+1),:],Y[:(i+1),0],k,jit=S[0,0])
        l2=L[(1.1,0.2,0.3+(i+1)*1e-8)].llk(i+1)[0]
        t2=time.clock()
        p= gpbo.core.GPdc.GP_LKonly(X[:(i+1),:],Y[:(i+1)],S[:(i+1)],D[:(i+1)],k).llk()
        #print(p,p2)
        t3=time.clock()
        assert(np.allclose(l1,p))
        #assert(np.allclose(l2,p))
        T.append([t1-t0,t2-t1,t3-t2])

    f,a = plt.subplots(1)
    n=len(T)
    a.plot(np.arange(1,n+1),[t[0] for t in T],'b')
    a.plot(np.arange(1,n+1),np.cumsum([t[0] for t in T]),'b--')
    a.plot(np.arange(1,n+1),[t[1] for t in T],'r')
    a.plot(np.arange(1,n+1),[t[2] for t in T],'g')
    a.set_yscale('log')
    a.set_xscale('log')
    f.savefig('figs/sctime.png')

ki = gpbo.core.GPdc.MAT52
Dim=2
k=gpbo.core.GPdc.kernel(ki,Dim,np.array([1.1,0.2,0.3]))
if __name__=="__main__":
    main()
