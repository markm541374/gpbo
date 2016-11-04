# Gaussian process with derivatives in c
#LK_only can only be used for log likelihood or posterior log likelihood
#GPcore takes a single hyperparameter value or a set. With a set the infer_ methods produce individual inferences while infer_post methods return the posterior
#cython: profile=True
#

__author__ = "mark"
__date__ = "$22-Nov-2015 21:27:19$"

#TODO draw doesn't work correclty on a GP created with no data
import numpy as np
cimport numpy as np
import scipy as sp
import ctypes as ct
import os
import sys
from copy import deepcopy as dc
from scipy.stats import norm as norms
from libc.math cimport log10, log, isnan
#print os.path.join(os.path.split(__file__)[0],'../../dist/Release/GNU-Linux/libGPshared.so')
from . import __file__ as fl
from scipy import linalg as spl

libGP = ct.cdll.LoadLibrary(os.path.join(os.path.split(fl)[0],'../cproj/libcproj.so')) #path to c-shared library
#print libGP
libGP.k.restype = ct.c_double

ctpd = ct.POINTER(ct.c_double)
cint = ct.c_int


class GP_LKonly:
    def __init__(self, X_s, Y_s, S_s, D_s, kf):
        cdef int n, D, i
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
    def __init__(self, X_s, Y_s, S_s, D_s, kf):
        #print [X_s, Y_s, S_s, D_s, kf]
        if type(kf) is kernel:
            self.size = 1
            kf = [kf]
        else:
            self.size = len(kf)
        allhyp = sp.hstack([k.hyp for k in kf])
        self.kf=kf
        [self.n ,self.D] = X_s.shape
        Dx = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_s]
        self.s = libGP.newGP_hypset(cint(self.D),cint(self.n),cint(kf[0].Kindex),X_s.ctypes.data_as(ctpd),Y_s.ctypes.data_as(ctpd),S_s.ctypes.data_as(ctpd),(cint*len(Dx))(*Dx),allhyp.ctypes.data_as(ctpd),cint(self.size))
        self.Y_s=Y_s
        
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_s]
        #print self.get_cho()
        libGP.presolv(self.s,cint(self.size))
        
        return
    
    def __del__(self):
        #print "here"
        libGP.killGP(self.s,cint(self.size))
        #print "also here"
        return
    
    def printc(self):
        print self.size
        libGP.ping(self.s, cint(self.size))
        return
    def get_cho(self):
        C=sp.empty([self.n,self.n*self.size])
        libGP.get_cho(self.s,cint(self.size),C.ctypes.data_as(ctpd))
        return C
    def infer_m(self,X_i,D_i):
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.vstack([sp.empty(ns)]*self.size)
        libGP.infer_m(self.s, cint(self.size), ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D),R.ctypes.data_as(ctpd))
        return R
    
    def infer_m_partial(self,X_i,D_i,ki,hyp):
        
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.vstack([sp.empty(ns)]*1)
        
        #libGP.infer_m_partial(self.s, cint(ki),hyp.ctypes.data_as(ctpd),ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D),R.ctypes.data_as(ctpd))
        libGP.infer_m_partial(self.s,cint(ki),hyp.ctypes.data_as(ctpd),ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D),R.ctypes.data_as(ctpd))
            
        
        return R
    
    def infer_m_post(self,X_i,D_i):
        ns=X_i.shape[0]
        R = self.infer_m(X_i,D_i)
        
        return sp.mean(R,axis=0).reshape([1,ns])
    
    
    def infer_full(self,X_i,D_i):
        cdef int ns,j,i
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**j for j in x])) for x in D_i]
        R=sp.vstack([sp.empty([ns+1,ns])]*self.size)
        libGP.infer_full(self.s, cint(self.size), ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D),R.ctypes.data_as(ctpd))
        m = sp.vstack([R[i*(ns+1),:] for i in range(self.size)])
        V = sp.vstack([R[(ns+1)*i+1:(ns+1)*(i+1),:] for i in range(self.size)])
        return [m,V]
    
    def infer_full_post(self,X_i,D_i):
        class MJMError(Exception):
            pass
        [m,V] = self.infer_full(X_i,D_i)
        cdef int i,ns
        ns=X_i.shape[0]
        cv = sp.zeros([ns,ns])

        for i in range(self.size):
            cv+=V[ns*i:ns*(i+1),:]
        cv= cv/self.size
        if self.size>1:
            cv+=sp.cov(m,rowvar=0,bias=1)
        #print "_________________"
        #print self.size
        #print sp.cov(m,rowvar=0,bias=1)
        #print V
        #print cv
        return [sp.mean(m,axis=0).reshape([1,ns]),cv]
    
    def infer_diag(self,X_i,D_i):
        cdef int i,j,ns
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**j for j in x])) for x in D_i]
        R=sp.vstack([sp.empty([2,ns])]*self.size)
        libGP.infer_diag(self.s,cint(self.size), ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D),R.ctypes.data_as(ctpd))

        m = sp.vstack([R[i*2,:] for i in range(self.size)])
        V = sp.vstack([R[i*2+1,:] for i in range(self.size)])
        return [m,V]
    
    def infer_diag_post(self,X_ii,D_i):
        cdef int ns
        X_i = dc(X_ii)
        ns = len(D_i)
        
        X_i.resize([ns,self.D])
        [m,V] = self.infer_diag(X_i,D_i)
        if sp.amin(V)<=-0.:
            class MJMError(Exception):
                pass
            print "negative/eq variance"
            print [m,V,X_i,D_i]
            print "_______________"
            #self.printc()
            raise(MJMError)
        if sp.amin(sp.var(m,axis=0))<-0.:
            class MJMError(Exception):
                pass
            print "negativevar of mean"
            print X_i.shape
            print [m,V,sp.var(m,axis=0),X_i,D_i]
            print "_______________"
            #self.printc()
            raise(MJMError)
        
        return [sp.mean(m,axis=0).reshape([1,ns]),(sp.mean(V,axis=0)+sp.var(m,axis=0)).reshape([1,ns])]
        
    
    def draw(self,X_i,D_i,m):
        #make m draws at X_i Nd, X, D, R, m
        cdef int i,ns
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([m*self.size,ns])
        libGP.draw(self.s, cint(self.size), cint(ns), X_i.ctypes.data_as(ctpd), (cint*len(D))(*D),R.ctypes.data_as(ctpd),cint(m))
        return R
    
    def draw_post(self,X_i,D_i,z):
        cdef int ns
        ns = X_i.shape[0]
        [m,V] = self.infer_full_post(X_i,D_i)
        R = sp.empty([ns,z])
        libGP.drawk(V.ctypes.data_as(ctpd),cint(ns),R.ctypes.data_as(ctpd),cint(z))
        R+=sp.hstack([m.T]*z)
     #R=sp.random.multivariate_normal(m.flatten(),V,z)
        return R.T
    
    def llk(self):
        R = sp.empty(self.size)
        libGP.llk(self.s, cint(self.size), R.ctypes.data_as(ctpd))
        return R
    
    def infer_LCB(self,X_i,D_i, p):
        cdef int ns,i
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([self.size,ns])
        libGP.infer_LCB(self.s, cint(self.size), ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D), ct.c_double(p), R.ctypes.data_as(ctpd))
        
        return R
    
    def infer_LCB_post(self,X_i,D_i,p):
        [m,v] = self.infer_diag_post(X_i,D_i)
        if sp.amin(v)<0.:
            print "negateive vriance: "
            print [m,v,X_i]
            #self.printc()
            class MJMError(Exception):
                pass
            raise MJMError()
        return m-p*sp.sqrt(v)
    
    def infer_EI(self,X_i,D_i):
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([self.size,ns])
        libGP.infer_EI(self.s, cint(self.size),ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D), R.ctypes.data_as(ctpd))
        return R
    
    def infer_EI_post(self,X_i,D_i):
        [m,v] = self.infer_diag_post(X_i,D_i)
        cdef int ns=X_i.shape[0]
        R=sp.empty([1,ns])
        cdef int i
        for i in range(ns):
            R[0,i] = EI(sp.amin(self.Y_s),m[0,i],v[0,i])
        
        return R
    
    def infer_lEI(self,X_i,D_i):
        ns=X_i.shape[0]
        D = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.empty([self.size,ns])
        libGP.infer_lEI(self.s, cint(self.size),ns,X_i.ctypes.data_as(ctpd),(cint*len(D))(*D), R.ctypes.data_as(ctpd))
        return R
#kf = gen_sqexp_k_d([1.,0.3])


SQUEXP = 0
LIN1 = 1
LINXSQUEXP = 2
LINSQUEXPXSQUEXP = 3
SQUEXP1SSQUEXP = 4
SSPS = 5
SQUEXPCS = 6
SQUEXPPS = 7
SQUEXPBS = 8
MAT52 = 9
MAT52CS = 10
MAT52PER = 11
MAT52PPT = 12
DEV=13
MATPP=14
class kernel(object):
    def __init__(self,K,D,H):
        self.dim = D
        self.hyp = sp.array(H)
        self.Kindex = K
        #ihyp are derived from the hyperparameters for speed and will be 1/h^2 etc.
        self.ihyp = sp.empty(self.hyp.shape[0])
        libGP.hypconvert(self.hyp.ctypes.data_as(ctpd),cint(self.Kindex), cint(self.dim), self.ihyp.ctypes.data_as(ctpd))
        return
    
    def __call__(self,x1, x2, d1=[sp.NaN], d2=[sp.NaN],gets=False):
        D1 = 0 if isnan(d1[0]) else int(sum([8**x for x in d1]))
        D2 = 0 if isnan(d2[0]) else int(sum([8**x for x in d2]))
        self.smodel=sp.empty(1)
        r=libGP.k(x1.ctypes.data_as(ctpd),x2.ctypes.data_as(ctpd), cint(D1),cint(D2),cint(self.dim),self.ihyp.ctypes.data_as(ctpd),cint(self.Kindex),self.smodel.ctypes.data_as(ctpd))
        if gets:
            return [r,self.smodel[0]]
        return r
    

class gen_sqexp_k_d():
    def __init__(self,theta):
        self.dim = len(theta)-1
        self.hyp = sp.array(theta)
        self.hypinv = sp.array([1./x**2 for x in theta])
        self.hypinv[0] = theta[0]**2
        self.Kindex = 0;
        return
    def __call__(self,x1, x2, d1=[sp.NaN], d2=[sp.NaN]):
        D1 = 0 if isnan(d1[0]) else int(sum([8**x for x in d1]))
        D2 = 0 if isnan(d2[0]) else int(sum([8**x for x in d2]))
        self.smodel=0.
        r=libGP.k(x1.ctypes.data_as(ctpd),x2.ctypes.data_as(ctpd), cint(D1),cint(D2),cint(self.dim),self.hypinv.ctypes.data_as(ctpd),cint(0),ct.byref(ct.c_double(self.smodel)))
        return r
    
class gen_lin1_k_d():
    def __init__(self,theta):
        self.hyp = sp.array(theta)
        self.Kindex = 1
        self.hypinv = sp.array(theta)
        self.hypinv[0] = self.hypinv[0]**2
        self.hypinv[2] = self.hypinv[2]**2
        
        return
    
    def __call__(self,x1, x2, d1=[sp.NaN], d2=[sp.NaN]):
        D1 = 0 if isnan(d1[0]) else int(sum([8**x for x in d1]))
        D2 = 0 if isnan(d2[0]) else int(sum([8**x for x in d2]))
        smodel = 0.
        r=libGP.k(x1.ctypes.data_as(ctpd),x2.ctypes.data_as(ctpd), cint(D1),cint(D2),cint(-42),self.hypinv.ctypes.data_as(ctpd),cint(1),ct.byref(ct.c_double(smodel)))
        return r

def searchMLEhyp(X,Y,S,D,lb,ub, ki, mx=20000,fg=-1e9):
    libGP.SetHypSearchPara(cint(mx),ct.c_double(fg))
    ns=X.shape[0]
    dim = X.shape[1]
    Dx = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D]
   
    hy = sp.empty(libGP.numhyp(cint(ki),cint(dim)))
    
    lk = sp.empty(1)
    r = libGP.HypSearchMLE(cint(dim),cint(len(Dx)),X.ctypes.data_as(ctpd),Y.ctypes.data_as(ctpd),S.ctypes.data_as(ctpd),(cint*len(Dx))(*Dx),lb.ctypes.data_as(ctpd),ub.ctypes.data_as(ctpd),cint(ki), hy.ctypes.data_as(ctpd),lk.ctypes.data_as(ctpd))
    
    return hy


def searchMAPhyp(X,Y,S,D,m,s, ki, MAPmargin = 1.8, mx=20000,fg=-1e9):
    libGP.SetHypSearchPara(cint(mx),ct.c_double(fg))
    ns=X.shape[0]
    dim = X.shape[1]
    Dx = [0 if isnan(x[0]) else int(sum([8**i for i in x])) for x in D]
    hy = sp.empty(libGP.numhyp(cint(ki),cint(dim)))
    
    lk = sp.empty(1)
    print "datasetsize = "+str(ns)
    r = libGP.HypSearchMAP(cint(dim),cint(len(Dx)),X.ctypes.data_as(ctpd),Y.ctypes.data_as(ctpd),S.ctypes.data_as(ctpd),(cint*len(Dx))(*Dx),m.ctypes.data_as(ctpd),s.ctypes.data_as(ctpd),ct.c_double(MAPmargin),cint(ki), hy.ctypes.data_as(ctpd),lk.ctypes.data_as(ctpd))
    #print "yyy"
    return hy

#just for quickly making test draws
def buildKsym_d(kf,x,d):
        #x should be  column vector
        (l,_)=x.shape
        K=sp.matrix(sp.empty([l,l]))
        cdef int i,j
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
