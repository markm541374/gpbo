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
lb = sp.array([0.]+[-1.]*(d-1))
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt,d,lb,ub,GPdc.SQUEXP,sp.array([1.5,0.55,0.25]))

kindex = GPdc.SQUEXP
mprior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
axis=0
value=0.
pesobj = PES.PES_inplane(X,Y,S,D,lb,ub,kindex,mprior,sprior,axis,value,DH_SAMPLES=8,DM_SAMPLES=8, DM_SUPPORT=400,DM_SLICELCBPARA=1.,AM_POLICY=PES.NOMIN,mode=ESutils.SUPPORT_SLICEEI)


def cfn(x):
    return 1.-(0.6*x[0])**0.1
def sfn(x):
    return 1e-4

[xmin,ymin,ierror] = pesobj.search_acq(cfn,sfn)
print [xmin,ymin,ierror]


np=151
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Sp = sp.array([[1e-3]*np]).T
Xp0 = sp.vstack([sp.array([i,xmin[1]]) for i in sup])
Xp1 = sp.vstack([sp.array([xmin[0],i]) for i in sup])

f,a = plt.subplots(d)
h0 = pesobj.query_pes(Xp0,Sp,Dp)
h1 = pesobj.query_pes(Xp1,Sp,Dp)
a[0].plot(sup,h0.flatten(),'b')
a[1].plot(sup,h1.flatten(),'b')

a0 = pesobj.query_acq(Xp0,Sp,Dp,cfn)
a1 = pesobj.query_acq(Xp1,Sp,Dp,cfn)
a[0].plot(sup,a0.flatten(),'r')
a[1].plot(sup,a1.flatten(),'r')

plt.show()