#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import ESutils
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
import GPdc
import PES

#-------------------------------------------------------------------------
#2d
nt=30
d=2
lb = sp.array([-1.]*d)
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt, d, lb, ub, GPdc.SQUEXP, sp.array([1.5, 0.35, 0.30]))

kindex = GPdc.SQUEXP
mprior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))

pesobj = PES.PES(X,Y,S,D,lb,ub,kindex,mprior,sprior,DH_SAMPLES=8,DM_SAMPLES=8, DM_SUPPORT=400,DM_SLICELCBPARA=1., mode=ESutils.SUPPORT_SLICEEI)

np=150
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Sp = sp.array([[1e-3]*np]).T
Xp0 = sp.vstack([sp.array([i,pesobj.Z[0,1]]) for i in sup])
Xp1 = sp.vstack([sp.array([pesobj.Z[0,0],i]) for i in sup])

f,a = plt.subplots(d)
h0 = pesobj.query_pes(Xp0,Sp,Dp)
h1 = pesobj.query_pes(Xp1,Sp,Dp)
a[0].plot(sup,h0.flatten(),'b')
a[1].plot(sup,h1.flatten(),'b')

def cfn(x,s):
    return ((x[0])**2)+x[1]+1.

a0 = pesobj.query_acq(Xp0,Sp,Dp,cfn)
a1 = pesobj.query_acq(Xp1,Sp,Dp,cfn)
a[0].plot(sup,a0.flatten(),'r')
a[1].plot(sup,a1.flatten(),'r')

cfn = lambda x,s:1.#s**-0.1
[xmin,ymin,ierror] = pesobj.search_acq(cfn,-4.,-1.)
print [xmin,ymin,ierror]

[xmin,ymin,ierror] = pesobj.search_pes(1e-4)
print [xmin,ymin,ierror]
plt.show()