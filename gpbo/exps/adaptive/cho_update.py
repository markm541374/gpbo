import numpy as np
import scipy as sp
from scipy.linalg import cho_factor
from scipy.linalg import cho_solve
import time
import tqdm
from matplotlib import pyplot as plt
import gpbo


class kupdate(object):
    def __init__(self,maxsize,x,y,k,jit=1e-6):
        self.k = k
        self.jit=jit
        self.C = np.zeros([maxsize,maxsize])
        self.X = np.zeros([maxsize,x.shape[1]])
        self.Y = np.zeros(maxsize)
        self.count=0
        for i in range(x.shape[0]):
            self.update(x[i,:],y[i])
        return

    def update(self,x,y):
        self.X[self.count,:]=x
        self.Y[self.count]=y
        for i in range(self.count):
            self.C[self.count,i]=(self.k(x,self.X[i,:])-self.C[self.count,:i].dot(self.C[i,:i]))/self.C[i,i]

        self.C[self.count,self.count] = np.sqrt(self.k(x,x)+self.jit-np.sum(self.C[self.count,:self.count]**2))

        self.count+=1

        #self.H = cho_factor(self.M[:self.count,:self.count],lower=True)
        #print(np.allclose(self.C[:self.count,:self.count],np.array(self.H[0])))
        #self.llk()
        return
    def llk(self):
        ym=np.mean(self.Y[:self.count])

        t2 = -0.5*(self.Y[:self.count]-ym).dot(cho_solve((self.C[:self.count,:self.count],True),self.Y[:self.count]-ym))
        t0 = -0.5*self.count*np.log(2*sp.pi)
        t1 = -np.sum(np.log(np.diag(self.C[:self.count,:self.count])))
        lk=t0+t1+t2

        #p= gpbo.core.GPdc.GP_LKonly(self.X[:self.count,:],self.Y[:self.count],self.jit*np.ones([self.count,1]),[[np.NaN]]*self.count,k).llk()
        #print(np.allclose(p,lk))
        return t0+t1+t2
def main():
    ki = gpbo.core.GPdc.MAT52
    Dim=2
    k=gpbo.core.GPdc.kernel(ki,Dim,np.array([1.1,0.2,0.3]))
    X = np.random.uniform(-1,1,[400,Dim])
    Y = sp.vstack([gpbo.core.objectives.shiftbraninojf(X[i,:],**{})[0] for i in range(X.shape[0])])
    S = np.ones_like(Y)*1e-6
    D = [[sp.NaN]]*Y.size

    G = kupdate(400,X[:0,:],Y[:0,0],k,jit=S[0,0])
    T=[]
    for i in tqdm.tqdm(range(100)):
        t0=time.clock()
        G.update(X[i,:],Y[i,0])
        l1=G.llk()
        t1=time.clock()
        Gf = kupdate(i+1,X[:(i+1),:],Y[:(i+1),0],k,jit=S[0,0])
        l2=Gf.llk()
        t2=time.clock()
        p= gpbo.core.GPdc.GP_LKonly(X[:(i+1),:],Y[:(i+1)],S[:(i+1)],D[:(i+1)],k).llk()
        t3=time.clock()
        assert(np.allclose(l1,l2))
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