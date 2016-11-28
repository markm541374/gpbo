# Gaussian process for 1d grid
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
import os
import sys
from copy import deepcopy as dc
from scipy.stats import norm as norms
from libc.math cimport log10, log, isnan
#print os.path.join(os.path.split(__file__)[0],'../../dist/Release/GNU-Linux/libGPshared.so')

from GPdc import kernel
class GP_LKonly:
    def __init__(self, X_s, Y_s, S_s, D_s, kf):
        raise NotImplementedError
        return
    
    def llk(self):

        raise NotImplementedError
        #log likelihood of data given hyperparametrs
        return self.l
    
    def plk(self,lm,ls):
        #log likelihood given lognormal prior over hyperparameters

        raise NotImplementedError
        return self.l+tmp

class GPcore:
    def __init__(self, X_s, Y_s, S_s, D_s, kf):
        assert type(kf) is kernel
        self.kf=kf
        self.n,d = X_s.shape
        assert(d==1)
        assert D_s==[[sp.NaN]]*self.n
        self.s = S_s[0,0]

        return
    
    def printc(self):
        raise NotImplementedError
        return

    def infer_m(self,X_i,D_i):
        raise NotImplementedError
        return R

    def infer_m_post(self,X_i,D_i):
        raise NotImplementedError
        return sp.mean(R,axis=0).reshape([1,ns])
    
    
    def infer_full(self,X_i,D_i):
        raise NotImplementedError
        return [m,V]
    
    def infer_full_post(self,X_i,D_i):
        raise NotImplementedError
        return [sp.mean(m,axis=0).reshape([1,ns]),cv]
    
    def infer_diag(self,X_i,D_i):
        raise NotImplementedError
        return [m,V]
    
    def infer_diag_post(self,X_ii,D_i):
        raise NotImplementedError
        return [sp.mean(m,axis=0).reshape([1,ns]),(sp.mean(V,axis=0)+sp.var(m,axis=0)).reshape([1,ns])]
        
    
    def draw(self,X_i,D_i,m):
        raise NotImplementedError
        return R
    
    def draw_post(self,X_i,D_i,z):
        raise NotImplementedError
        return R.T
    
    def llk(self):
        raise NotImplementedError
        return R
    
    def infer_LCB(self,X_i,D_i, p):
        raise NotImplementedError

        return R
    
    def infer_LCB_post(self,X_i,D_i,p):
        raise NotImplementedError
        return m-p*sp.sqrt(v)
    
    def infer_EI(self,X_i,D_i):
        raise NotImplementedError
        return R
    
    def infer_EI_post(self,X_i,D_i,wrt=False):
        raise NotImplementedError
        return R
    
    def infer_lEI(self,X_i,D_i):
        raise NotImplementedError
        return R
#kf = gen_sqexp_k_d([1.,0.3])



def searchMLEhyp(X,Y,S,D,lb,ub, ki, mx=20000,fg=-1e9):
    raise NotImplementedError
    return hy


def searchMAPhyp(X,Y,S,D,m,s, ki, MAPmargin = 1.8, mx=20000,fg=-1e9):
    raise NotImplementedError
    return hy


