#!/usr/bin/env python2
#encoding: UTF-8

# Comparing var fid PES to fixed PES and EI on branin for various fixed step sizes. cost/noise fixed at (1e-6)/s)**0.2


import scipy as sp
from matplotlib import pyplot as plt
import gpbo.core.GPdc
import gpbo.core.optutils
import gpbo.core.search
import os
import pickle
import gpbo.core.wrappingLogistic
import time

print 'start'

d=4
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])


def cfn(s):
    return 1.

def ojf2(x,s,d,override=False):
    x.resize([1,x.size])
    f = (0.1-x[0,0])**2+(0.2-x[0,1])**2+(-0.1-x[0,2])**2+(-0.2-x[0,3])**2
    return [f,1.]

truelb = [0.,0.,20,5]
trueub = [10,1,2000,2000]

def ojf(x,s,d,override=False):
    x.resize([1,x.size])
    xscale = [0.,0.,0.,0.]
    for i in xrange(4):
        xscale[i] = truelb[i]+(x[0,i]+1.)*0.5*(trueub[i]-truelb[i])
    t0 = time.clock()
    f = gpbo.core.wrappingLogistic.main({'lrate':xscale[0], 'l2_reg':xscale[1], 'batchsize':xscale[2], 'n_epochs':80}, fold=1, folds=1, downsize=1)
    t1 = time.clock()
    return [f,t1-t0]


braninmin = 0.
kindex = gpbo.core.GPdc.MAT52
prior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
kernel = [kindex,prior,sprior]
nreps = 2
bd = 4*60
slist = [1e-9]
print 'start'
f,a = plt.subplots(3)
import os

for s in slist:
    
    names = ["../cache/mnist/EIMLE_5mnist"+str(int(100*sp.log10(s)))+"_"+"_"+str(i)+".p" for i in xrange(nreps)]
    results = gpbo.core.search.multiMLEFS(ojf, lb, ub, kernel, s, bd, names)
    yr = [r[11].flatten() for r in results]
    C = [r[5] for r in results]
    
for i in xrange(nreps):
    m = yr[i]
    a[2].plot([sum(C[i][:j]) for j in xrange(len(C[i]))],m.flatten(),'x-')
    a[1].plot(C[i],'r')


print "reccomend: "+str([r[4][-1,:] for r in results])

f.savefig("../figs/braninnoise.png")
plt.show()