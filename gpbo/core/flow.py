# Gaussian process with derivatives in c
#LK_only can only be used for log likelihood or posterior log likelihood
#GPcore takes a single hyperparameter value or a set. With a set the infer_ methods produce individual inferences while infer_post methods return the posterior
#cython: profile=True
#
from __future__ import print_function
xrange=range
__author__ = "mark"
__date__ = "$22-Nov-2015 21:27:19$"

import numpy as np
import scipy as sp
import copy
import os
from scipy.stats import norm as norms
from numpy import log10, log, isnan, exp
from gpbo.core import flowkern
import GPflow as gpf

class flow_Error(Exception):
    pass

class GP_LKonly:
    def __init__(self, X_s, Y_s, S_s, D_s, kf):
        raise NotImplementedError
        n ,D = X_s.shape
        self.hyp = kf.hyp
        R = ct.c_double()
        Dc = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_s]
        libGP.newGP_LKonly(cint(D),cint(n),X_s.ctypes.data_as(ctpd),Y_s.ctypes.data_as(ctpd),S_s.ctypes.data_as(ctpd),(cint*len(Dc))(*Dc), cint(kf.Kindex), kf.hyp.ctypes.data_as(ctpd),ct.byref(R))
        self.l = R.value
        
        return
    
    def llk(self):
        #log likelihood of data given hyperparametrs
        return self.l
    
    def plk(self,lm,ls):
        #log likelihood given lognormal prior over hyperparameters
        tmp = 0.

        for i,h in enumerate(self.hyp):
            tmp -= 0.5*((log10(h)-lm[i])**2)/ls[i]**2
        
        return self.l+tmp

class GPcore:
    def __init__(self, X, Y, S, D, kf):
        self.X=X
        self.Y=Y
        if isinstance(kf,kernel):
            self.d = kf.dim
            self.size=1
            kf = [kf]
        else:
            self.size=len(kf)
            self.d = kf[0].dim
        x = np.hstack([X,np.array([[0 if isnan(x[0]) else (sum([1. if i==j else 0 for i in x ])) for x in D] for j in range(self.d)]).T])
        X = np.hstack([x,S.reshape(x.shape[0],1)])
        #self.m = [gpf.gpr.GPR(X,Y,klist[kf[0].Kindex](kf[0].dim)) for i in range(self.size)]
        self.m = [gpf.gpr.GPR(X,Y,kf[i].K) for i in range(self.size)]
        for i in range(self.size):
            if kf[i].hyparray:
                self.m[i].likelihood.variance.transform = gpf.transforms.Log1pe(lower=0.001*np.min(S))
                self.m[i].likelihood.variance=0.0011*np.min(S)
                self.m[i].likelihood.variance.fixed= True
            else:
                self.m[i].set_parameter_dict(kf[i].hyp)
        return

    def printc(self):
        for i in range(self.size):
            print(self.m[i])
        return

    def infer_m(self,X, D):
        x = np.hstack([X,np.array([[0 if isnan(x[0]) else (sum([1. if i==j else 0 for i in x ])) for x in D] for j in range(self.d)]).T])
        m = np.empty([X.shape[0],self.size])
        for i in range(self.size):
            m[:,i:i+1], _ = self.m[i].predict_f(np.hstack([x,np.zeros(shape=x.shape)]))
        return m.T

    def infer_m_post(self,X_,D_i):
        X_i = copy.copy(X_)
        ns=X_i.shape[0]
        R = self.infer_m(X_i,D_i)
        return sp.mean(R,axis=0).reshape([1,ns])

    def infer_full(self,X,D):
        x = np.hstack([X,np.array([[0 if isnan(x[0]) else (sum([1. if i==j else 0 for i in x ])) for x in D] for j in range(self.d)]).T])
        n = X.shape[0]
        m = np.empty([n,self.size])
        V = np.empty([n,n*self.size])
        for i in range(self.size):
            m_ ,V_ = self.m[i].predict_f_full_cov(np.hstack([x,np.zeros(shape=x.shape)]))
            m[:,i:i+1] = m_
            V[:,i*n:(i+1)*n] = np.squeeze(V_)
        return m.T,V.T

    def infer_full_post(self,X_,D_i):
        X_i = copy.copy(X_)
        [m,V] = self.infer_full(X_i,D_i)
        ns=X_i.shape[0]
        cv = sp.zeros([ns,ns])

        for i in range(self.size):
            cv+=V[ns*i:ns*(i+1),:]
        cv= cv/self.size
        if self.size>1:
            cv+=sp.cov(m,rowvar=0,bias=1)
        return [sp.mean(m,axis=0).reshape([1,ns]),cv]

    def infer_diag(self, X, D):
        x = np.hstack([X,np.array([[0 if isnan(x[0]) else (sum([1. if i==j else 0 for i in x ])) for x in D] for j in range(self.d)]).T])

        m = np.empty([X.shape[0],self.size])
        V = np.empty([X.shape[0],self.size])
        for i in range(self.size):
            m[:,i:i+1], V[:,i:i+1] = self.m[i].predict_f(np.hstack([x,np.zeros(shape=x.shape)]))
        if sp.amin(V)<=-0.:
            print( "negative/eq variance")
            print( "_______________")
            print(V)
            raise(flow_Error)
        return m.T, V.T

    def infer_diag_post(self,X_,D_i):
        X_i = copy.copy(X_)
        ns = len(D_i)

        X_i.resize([ns,self.d])
        [m,V] = self.infer_diag(X_i,D_i)
        if sp.amin(V)<=-0.:
            print( "negative/eq variance")
            print( "_______________")
            raise(flow_Error)
        if sp.amin(sp.var(m,axis=0))<-0.:
            print( "negativevar of mean")
            print( "_______________")
            raise(flow_Error)

        return [sp.mean(m,axis=0).reshape([1,ns]),(sp.mean(V,axis=0)+sp.var(m,axis=0)).reshape([1,ns])]

    def draw(self,X,D,n):
        x = np.hstack([X,np.array([[0 if isnan(x[0]) else (sum([1. if i==j else 0 for i in x ])) for x in D] for j in range(self.d)]).T])
        s = np.empty([self.size*n,X.shape[0]])
        for i in range(self.size):
            s[n*i:n*(i+1),:] = np.squeeze(self.m[i].predict_f_samples(np.hstack([x,np.zeros(shape=x.shape)]),n))
        return s

    def draw_post(self,X,D,n):
        x = np.hstack([X,np.array([[0 if isnan(x[0]) else (sum([1. if i==j else 0 for i in x ])) for x in D] for j in range(self.d)]).T])
        s = np.empty([n,X.shape[0]])
        J = np.random.randint(0,self.size,n)
        for i in range(n):
            try:
                s[i,:]= np.squeeze(self.m[J[i]].predict_f_samples(np.hstack([x,np.zeros(shape=x.shape)]),1))
            except:
                K= np.squeeze(self.m[J[i]].predict_f_full_cov(np.hstack([x,np.zeros(shape=x.shape)]))[1])
                print(K.shape)
                print(K)
                from scipy.linalg import cho_factor
                for i in range(-12,0):
                    try:
                        cho_factor(K+sp.eye(K.shape[0])*10**i)
                        print(i)
                    except:
                        pass
                raise
        return s

    def llk(self):
        return np.array([self.m[i].compute_log_likelihood() for i in range(self.size)])

    def infer_LCB(self,X, D, p):
        m,v = self.infer_diag(X, D)
        return m-p*np.sqrt(v)

    def infer_LCB_post(self,X, D, p):
        m,v = self.infer_diag_post(X, D)
        return m-p*np.sqrt(v)
        raise NotImplementedError
        X_i = copy.copy(X_)
        [m,v] = self.infer_diag_post(X_i,D_i)
        if sp.amin(v)<0.:
            print( "negateive vriance: ")
            #print( [m,v,X_i])
            #self.printc()
            class GPdcError(Exception):
                pass
            raise GPdcError()
        return m-p*sp.sqrt(v)
    
    def infer_EI(self,X,D,fixI=False,I=0.):
        m,v = self.infer_diag(X,D)
        if not fixI:
            I=np.infty
            for i in range(len(self.Y)):
                if sum(X[i,self.d:])==0:
                    I = min(I,self.Y[i,0])
        E = sp.empty([self.size,X.shape[0]])
        for i in range(self.size):
            E[i,:] = EI(m[i,:],sp.sqrt(v[i,:]),I)
        return E

    def infer_EI_post(self,X,D,fixI=False,I=0.):
        E = self.infer_EI(X, D, fixI=fixI, I=I)
        ns = X.shape[0]
        return np.mean(E, axis=0).reshape([1, ns])

    def infer_lEI(self,X,D,fixI=False,I=0.):
        print('xxxx')
        m,v = self.infer_diag(X,D)
        if not fixI:
            I=np.infty
            for i in range(len(self.Y)):
                if sum(X[i,self.d:])==0:
                    I = min(I,self.Y[i,0])
        E = sp.empty([self.size,X.shape[0]])
        for i in range(self.size):
            E[i,:] = lEI(m[i,:],sp.sqrt(v[i,:]),I)
        return E

    def infer_lEI_post(self,X_,D_i,fixI=False,I=0.):
        E = self.infer_lEI(X_,D_i,fixI=fixI,I=I)
        ns=X_.shape[0]
        return sp.log(sp.nanmean(sp.exp(E),axis=0)).reshape([1,ns])




SQUEXP = 0
def squexphs(d,h):
    k0 = flowkern.dkern(d)
    k1 = flowkern.Pointwise_Hetroskedastic(1,active_dims=[flowkern.DXmax*d])
    k0.lengthscales=h[1:1+d]
    k0.variance=h[0]**2
    return k0+k1

klist = [squexphs, flowkern.dkern]
class kernel(object):
    def __init__(self,K,D,H):
        self.dim = D
        if not isinstance(H,dict):
            self.hyp = sp.array(H)
            self.hyparray=True
        else:
            self.hyparray=False
            self.hyp = H
        self.Kindex = K
        self.K = klist[K](D,H)
        return
    
    def __call__(self,x1, x2, d1=[sp.NaN], d2=[sp.NaN],gets=False):
        raise NotImplementedError
        return
    



def searchMLEhyp(X,Y,S,D,lb,ub, ki, mx=20000,fg=-1e9):
    g = GPcore(X,Y,S,D,kernel(ki,X.shape[1],np.ones(len(lb))))
    g.m[0].optimize()
    p = g.m[0].get_parameter_dict()
    print('MLE hyperparameters:\n{}'.format(p))
    return p


def searchMAPhyp(X,Y,S,D,m,s, ki, MAPmargin = 1.8, mx=20000,fg=-1e9):
    g = GPcore(X,Y,S,D,kernel(ki,X.shape[1],np.ones(len(m))))
    g.m[0].kern.lengthscales.prior = gpf.priors.LogNormal(np.log(10)*m[1:],(np.log(10)*s[1:])**2)
    g.m[0].kern.variance.prior = gpf.priors.LogNormal(np.log(10)*m[0],(np.log(10)*s[0])**2)
    g.m[0].optimize()
    p = g.m[0].get_parameter_dict()
    print('MAP hyperparameters:\n{}'.format(p))
    return p

#just for quickly making test draws
def buildKsym_d(kf,x,d):
        raise NotImplementedError
        #x should be  column vector
        (l,_)=x.shape
        K=sp.matrix(sp.empty([l,l]))
        for i in range(l):
            K[i,i]=kf(x[i,:],x[i,:],d1=d[i],d2=d[i])+10**-10
            for j in range(i+1,l):
                K[i,j]=kf(x[i,:],x[j,:],d1=d[i],d2=d[j])
                
                K[j,i]=K[i,j]
                
        return K
    
def EI(ER,mu,sigma):
        alpha=(-ER+mu)/sigma
        Z = norms.cdf(-alpha)
        
        if Z==0.0:
            return sp.matrix(0.0)
	#print "alpha: "+str(alpha)
	#print "Z: "+str(Z)
        E=-mu+norms.pdf(alpha)*sigma/Z+ER
        ei=Z*E
        if np.isfinite(ei):
            return sp.matrix(ei)
        else:
            return sp.matrix(0.0)

def draw(m_,V_,z):
        raise NotImplementedError
        ns = V_.shape[0]
        m = sp.array([[i for i in (m_)]])
        V = copy.copy(V_)
        R = sp.empty([ns,z])
        libGP.drawk(V.ctypes.data_as(ctpd),cint(ns),R.ctypes.data_as(ctpd),cint(z))
        R+=sp.hstack([m.T]*z)
     #R=sp.random.multivariate_normal(m.flatten(),V,z)
        return copy.copy(R).T
from scipy import stats
SQRT_1_2PI = 1/np.sqrt(2*np.pi)

def EI(m, s, y ):
    N = len(m)
    R = sp.empty(N)
    for i in range(N):
        S = (y-m[i])/s[i]
        c = stats.norm.cdf(S)
        p = stats.norm.pdf(S)
        R[i] = (y-m[i])*c+s[i]*p
    return R

def lEI(m, s, y ):
    N=len(m)
    R = EI(m,s,y)
    for i in range(N):
        #TODO this switch to log puts a kink in the curve
        if (R[i]<=0.):

            S = (y-m[i])/s[i]
            R[i] = np.log(s[i])+np.log(SQRT_1_2PI) - 0.5*S**2
            #print("logversion: %f\n",R[i]);
        else:
            #print("regversion: %f %f\n",R[i], log(R[i]));
            R[i] = np.log(R[i])
    return R

