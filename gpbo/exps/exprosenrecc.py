#!/usr/bin/env python2
#encoding: UTF-8

# for EI looking at the difference between the posterior min point and the


import scipy as sp
from matplotlib import pyplot as plt
import gpbo.core.GPdc
import gpbo.core.OPTutils
import gpbo.core.search
import os
import pickle


d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])
pwr = 0.2
cfn = lambda s:((1e-6)/s)**pwr
ojf = gpbo.core.OPTutils.genbanana(cfn=cfn)
rosenmin = 0.
kindex = gpbo.core.GPdc.MAT52
prior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
kernel = [kindex,prior,sprior]
nreps = 12
bd = 100
s = 1e-6
f,a = plt.subplots(2)

names = ["../cache/rosenrecc/EIMLE_"+str(int(100*sp.log10(s)))+"_"+str(pwr)+"_"+str(i)+".p" for i in xrange(nreps)]
results = gpbo.core.search.multiMLEFS(ojf, lb, ub, kernel, s, bd, names)

C = results[0][5]

yr = [r[11].flatten() for r in results]    
Z = sp.vstack(yr)-rosenmin
m = sp.mean(sp.log10(Z),axis=0)
v = sp.var(sp.log10(Z),axis=0)
sq = sp.sqrt(v)
a[1].fill_between(sp.array([sum(C[:j]) for j in xrange(len(C))]),(m-sq).flatten(),(m+sq).flatten(),facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
a[1].plot([sum(C[:j]) for j in xrange(len(C))],m.flatten(),'b.-')

yy = [r[10].flatten() for r in results]    
Zy = sp.vstack(yy)-rosenmin
my = sp.mean(sp.log10(Zy),axis=0)
vy = sp.var(sp.log10(Zy),axis=0)
sqy = sp.sqrt(vy)
a[1].fill_between(sp.array([sum(C[:j]) for j in xrange(len(C))]),(my-sqy).flatten(),(my+sqy).flatten(),facecolor='lightgreen',edgecolor='lightgreen',alpha=0.5)
a[1].plot([sum(C[:j]) for j in xrange(len(C))],my.flatten(),'g.-')
a[1].set_ylabel("log10 regret")


for i in xrange(nreps):
    a[0].plot([sum(C[:j]) for j in xrange(len(C))],sp.log10(Z[i,:]).flatten(),'b-')
    a[0].plot([sum(C[:j]) for j in xrange(len(C))],sp.log10(Zy[i,:]).flatten(),'g-')
f.savefig("../figs/rosenrecc.png")
plt.show()

