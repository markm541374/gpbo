#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from scipy import stats as sps
from scipy import linalg as spl
import scipy as sp
from matplotlib import pyplot as plt
import ESutils

import GPdc

nt=80
d=1
lb = sp.array([-1.]*d)
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt, d, lb, ub, GPdc.SQUEXP, sp.array([0.9, 0.25]), s=0.)
S*=0.
for i in xrange(nt):
    x = X[i,0]
    s = -(1e-2)*(x-1.)*(x+1.1)
    Y[i,0]+= sps.norm.rvs(0,sp.sqrt(s))

f0 = plt.figure()
a0 = plt.subplot(111)
a0.plot(sp.array(X[:,0]).flatten(),Y,'g.')

lb = sp.array([-2.,-2.,-9,-2.,-2.])
ub = sp.array([2.,2.,-1,2.,2.])
MLEH =  GPdc.searchMLEhyp(X, Y, S, D, lb, ub, GPdc.SQUEXPPS, mx=20000)

mprior = sp.array([0.,-1.,-5.,-0.5,0.5])
sprior = sp.array([1.,1.,3.,1.,1.])

MAPH = GPdc.searchMAPhyp(X, Y, S, D, mprior, sprior, GPdc.SQUEXPPS, mx=20000)
print "MLEH: "+str(MLEH)
print "MAPH: "+str(MAPH)
G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXPPS, 1, sp.array(MAPH)))



print G.llk()

np=180
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,v] = G.infer_diag(Xp,Dp)
sq = sp.sqrt(v)
S = sp.empty([np,1])
for i in xrange(np):
    S[i,0] = -MAPH[2]*(Xp[i,0]-MAPH[3])*(Xp[i,0]-MAPH[4])
sc= sp.sqrt(S.flatten())

a0.plot(sup,m.flatten())
sc = sp.sqrt(S.flatten())

a0.fill_between(sup, sp.array(m-1.*sq).flatten(), sp.array(m+1.*sq).flatten(), facecolor='lightblue',edgecolor='lightblue')
a0.plot(sup,(m+2*sc).flatten(),'g')
a0.plot(sup,(m-2*sc).flatten(),'g')
plt.show()
