#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import lkojf
import scipy as sp
import time
from matplotlib import pyplot as plt
import GPdc

L = lkojf.lkojf()
nf = 400
L.makedata_default(nf)
print L.llk(sp.array([1.3,0.13,0.2,1e-2]))#


np = 200
X = sp.array([sp.linspace(0,nf,np)]).T

H = sp.empty([np,1])
T = sp.empty([np,1])
for i in xrange(np):
    [H[i,0],T[i,0]] = L.llks(sp.array([1.3,0.13,0.2,1e-2]),int(X[i,0]))


#lb = sp.array([0.,0.,-4.,-0.5*float(nf),float(nf)])
#ub = sp.array([4.,3.,3.,0.,1.5*float(nf)])
#MLEH =  GPdc.searchMLEhyp(X,H,sp.zeros([np,1]),[[sp.NaN]]*(np),lb,ub,GPdc.SQUEXPPS,mx=10000)
#G = GPdc.GPcore(X.copy(),H,sp.zeros([np,1]),[[sp.NaN]]*(np),GPdc.kernel(GPdc.SQUEXPPS,1,sp.array(MLEH)))
lb = sp.array([0.,0.,-4.,-1.,-3.])
ub = sp.array([4.,3.,3.,0.5,3.])
MLEH =  GPdc.searchMLEhyp(1. / float(nf) * X, H, sp.zeros([np, 1]), [[sp.NaN]] * (np), lb, ub, GPdc.SQUEXPBS, mx=10000)
G = GPdc.GPcore(1. / float(nf) * X.copy(), H, sp.zeros([np, 1]), [[sp.NaN]] * (np), GPdc.kernel(GPdc.SQUEXPBS, 1, sp.array(MLEH)))

[m,v] = G.infer_diag(1./float(nf)*X,[[sp.NaN]]*(np))

S = sp.empty([np,1])
for i in xrange(np):
    S[i,0] = MLEH[2]*((1./float(nf)*X[i,0])**(MLEH[3]*MLEH[4]))* ((1.-1./float(nf)*X[i,0])**(MLEH[3]*(1.-MLEH[4])))

lbc = sp.array([-4.,0.,-6.])
ubc = sp.array([2.,3.,0.])
MLEC =  GPdc.searchMLEhyp(X, sp.log(T), sp.zeros([np, 1]), [[sp.NaN]] * (np), lbc, ubc, GPdc.SQUEXPCS, mx=10000)
C = GPdc.GPcore(X.copy(), sp.log(T), sp.zeros([np, 1]), [[sp.NaN]] * (np), GPdc.kernel(GPdc.SQUEXPCS, 1, sp.array(MLEC)))

f,a = plt.subplots(2)
a[0].plot(X.flatten(),H.flatten(),'g.')
s = sp.sqrt(v).flatten()
a[0].fill_between(X.flatten(),m.flatten()-2*s,m.flatten()+2*s, facecolor='lightblue',edgecolor='lightblue')
a[0].plot(X.flatten(),m.flatten(),'b')
a[0].plot(X.flatten(),m.flatten()+2*sp.sqrt(S.flatten()),'g')
a[0].plot(X.flatten(),m.flatten()-2*sp.sqrt(S.flatten()),'g')

a[1].plot(X.flatten(),T.flatten(),'g.')
[mc,vc] = C.infer_diag(X,[[sp.NaN]]*(np))
sc = sp.sqrt(vc).flatten()
a[1].fill_between(X.flatten(),sp.exp(mc.flatten()-2*sc),sp.exp(mc.flatten()+2*sc), facecolor='lightblue',edgecolor='lightblue')
a[1].plot(X.flatten(),sp.exp(mc.flatten()),'b')
plt.show()