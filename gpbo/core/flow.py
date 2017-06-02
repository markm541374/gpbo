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
        if isinstance(kf,kernel):
            self.d = kf.dim
            self.size=1
            kf = [kf]
        else:
            self.size=len(kf)
            self.d = kf[0].dim
        x = np.hstack([X,np.array([[0 if isnan(x[0]) else (sum([1. if i==j else 0 for i in x ])) for x in D] for j in range(self.d)]).T])
        self.m = [gpf.gpr.GPR(x,Y,klist[kf[0].Kindex](kf[0].dim)) for i in range(self.size)]
        for i in range(self.size):
            if kf[i].hyparray:
                self.m[i].kern.lengthscales=kf[i].hyp[1:1+self.d]
                self.m[i].kern.variance=kf[i].hyp[0]**2
                self.m[i].likelihood.variance=0.25
                self.m[i].likelihood.variance.fixed= True
            else:
                self.m[i].set_parameter_dict(kf[i].hyp)
        return

    def printc(self):
        for i in range(self.size):
            print(self.m[i])
        return

    #def get_cho(self):
     #   raise NotImplementedError
    #    C=sp.empty([self.n,self.n*self.size])
   #     libGP.get_cho(self.s,cint(self.size),C.ctypes.data_as(ctpd))
    #    return C

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

    def draw(self,X_i,D_i,m):
        raise NotImplementedError
        #make m draws at X_i Nd, X, D, R, m
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([m*self.size,ns])
        libGP.draw(self.s, cint(self.size), cint(ns), X_i.ctypes.data_as(ctpd), (cint*len(D))(*D),R.ctypes.data_as(ctpd),cint(m))
        return R
    
    def draw_post(self,X_i,D_i,z):
        raise NotImplementedError
        ns = X_i.shape[0]
        [m,V] = self.infer_full_post(X_i,D_i)
        R = sp.empty([ns,z])
        libGP.drawk(V.ctypes.data_as(ctpd),cint(ns),R.ctypes.data_as(ctpd),cint(z))
        R+=sp.hstack([m.T]*z)
     #R=sp.random.multivariate_normal(m.flatten(),V,z)
        return R.T
    
    def llk(self):
        return np.array([self.m[i].compute_log_likelihood() for i in range(self.size)])

    def infer_LCB(self,X_,D_i, p):
        raise NotImplementedError
        X_i = copy.copy(X_)
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([self.size,ns])
        libGP.infer_LCB(self.s, cint(self.size), ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D), ct.c_double(p), R.ctypes.data_as(ctpd))
        
        return R
    
    def infer_LCB_post(self,X_,D_i,p):
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
    
    def infer_EI(self,X_,D_i,fixI=False,I=0.):
        raise NotImplementedError
        X_i = copy.copy(X_)
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([self.size,ns])

        libGP.infer_EI(self.s, cint(self.size),ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D), R.ctypes.data_as(ctpd),cbool(fixI),cdbl(I))
        return R
    
    def infer_EI_post(self,X_,D_i,fixI=False,I=0.):
        raise NotImplementedError
        E = self.infer_EI(X_,D_i,fixI=fixI,I=I)
        ns=X_.shape[0]

        return sp.mean(E,axis=0).reshape([1,ns])

    def infer_lEI(self,X_,D_i,fixI=False,I=0.):
        raise NotImplementedError
        X_i = copy.copy(X_)
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([self.size,ns])
        libGP.infer_lEI(self.s, cint(self.size),ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D), R.ctypes.data_as(ctpd),cbool(fixI),cdbl(I))
        return R

    def infer_lEI_post(self,X_,D_i,fixI=False,I=0.):
        raise NotImplementedError
        E = self.infer_lEI(X_,D_i,fixI=fixI,I=I)
        ns=X_.shape[0]
        #print(E)
        #print(sp.log(sp.nanmean(sp.exp(E),axis=0)))
        return sp.log(sp.nanmean(sp.exp(E),axis=0)).reshape([1,ns])




SQUEXP = 0

klist = [flowkern.dkern]
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