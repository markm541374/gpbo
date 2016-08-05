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

nt=22
d=1
lb = sp.array([-1.]*d)
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt,d,lb,ub,GPdc.SQUEXP,sp.array([0.9,0.25]),s=1e-8)
S*=0.
f0 = plt.figure()
a0 = plt.subplot(111)
a0.plot(sp.array(X[:,0]).flatten(),Y,'g.')

lb = sp.array([-2.,-2.,-9])
ub = sp.array([2.,2.,-1])
MLEH =  GPdc.searchMLEhyp(X,Y,S,D,lb,ub,GPdc.SQUEXPCS,mx=10000)

mprior = sp.array([0.,-1.,-5.])
sprior = sp.array([1.,1.,3.])

MAPH = GPdc.searchMAPhyp(X,Y,S,D,mprior,sprior,GPdc.SQUEXPCS,mx=10000)
print "MLEH: "+str(MLEH)
print "MAPH: "+str(MAPH)
G = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXPCS,1,sp.array(MLEH)))


print G.llk()

np=180
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,v] = G.infer_diag(Xp,Dp)
a0.plot(sup,m.flatten())
sq = sp.sqrt(v)

a0.fill_between(sup, sp.array(m-1.*sq).flatten(), sp.array(m+1.*sq).flatten(), facecolor='lightblue',edgecolor='lightblue')

plt.show()
