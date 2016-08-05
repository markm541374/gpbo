#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
import GPdc
import ESutils
nt=2
d=1
lb = sp.array([-1.])
ub = sp.array([1.])
[X,Y,S,D] = ESutils.gen_dataset(nt,d,lb,ub,GPdc.SQUEXP,sp.array([1.5,0.35]),s=1e-3)

g = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,1,sp.array([1.5,0.35])))

np=100
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])
[m,V] = g.infer_diag_post(Xp,Dp)

f,a = plt.subplots(3)
s = sp.sqrt(V[0,:])
a[0].fill_between(sup,sp.array(m[0,:]-2.*s).flatten(),sp.array(m[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[0].plot(sup,m[0,:].flatten())
a[0].plot(sp.array(X[:,0]).flatten(),Y,'g.')



ns=500
R = ESutils.draw_support(g,lb,ub,ns,ESutils.SUPPORT_UNIFORM)
w = sp.ones(ns)
R2 = ESutils.draw_support(g,lb,ub,ns,ESutils.SUPPORT_SLICELCB,para=2.)
w2 = sp.exp(-g.infer_LCB_post(R2,[[sp.NaN]]*ns,2.))

a[0].twinx().plot(R2.flatten(),w2.flatten(),'rx')

[m0,V0] = g.infer_full_post(R,[[sp.NaN]]*ns)
ndr=10000
Z = g.draw_post(R, [[sp.NaN]]*ns,ndr)
Z2 = g.draw_post(R2, [[sp.NaN]]*ns,ndr)
H = sp.zeros(ns)
H2 = sp.zeros(ns)
for i in xrange(ndr):
    H[sp.argmin(Z[i,:])]+=1
    H2[sp.argmin(Z2[i,:])]+=1
    
a[1].plot(R,H.flatten(),'bx')
a[2].plot(R2,H2.flatten(),'rx')
a[1].twinx().plot(R2,H2.flatten()*w2.flatten(),'r.')
plt.show()