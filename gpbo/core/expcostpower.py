#!/usr/bin/env python2
#encoding: UTF-8

# Comparing how pes does on branin with various pwers to vf.


import scipy as sp
from matplotlib import pyplot as plt
import GPdc
import OPTutils
import search
import os
import pickle
print 'start'
print 'start'
d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])
pwr = 0.3
def cfn(s):
    #print s
    #print 'cfn'+str(((1e-6)/s)**pwr)
    return ((1e-6)/s)**pwr
ojf = OPTutils.genbranin(cfn=cfn)
braninmin = 0.39788735772973816
kindex = GPdc.MAT52
prior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
kernel = [kindex,prior,sprior]


nreps = 2
bd = 10
slist = [1e-4, 1e-5,1e-6,1e-7]
plist = [0.1,0.2,0.3,0.4,0.5,0.6]
print 'start'
lines = [[[],[],[],[]] for s in slist]
pline = [[],[],[],[]]
for l,pwr in enumerate(plist):
    f,a = plt.subplots(len(slist)+1)
    
    f.canvas.set_window_title(str(pwr))
    def cfn(s):
    
        print 'cfn'+str(((1e-6)/s)**pwr)
        return ((1e-6)/s)**pwr
    for k,s in enumerate(slist):
    
        names = ["../cache/costpower/PESFS_"+str(int(100*sp.log10(s)))+"_"+str(pwr)+"_"+str(i)+".p" for i in xrange(nreps)]
        results = search.multiPESFS(ojf,lb,ub,kernel,s,bd,names)
        zr = [r[11].flatten() for r in results]
        C = results[0][5]


        Z = sp.vstack(zr)-braninmin
        m = sp.mean(sp.log10(Z),axis=0)
        v = sp.var(sp.log10(Z),axis=0)
        sq = sp.sqrt(v)
        a[k].fill_between(sp.array([sum(C[:j]) for j in xrange(len(C))]),(m-sq).flatten(),(m+sq).flatten(),facecolor='lightgreen',edgecolor='lightgreen',alpha=0.5)
        a[k].plot([sum(C[:j]) for j in xrange(len(C))],m.flatten(),'.-')

        lines[k][0].append(pwr)
        lines[k][1].append(m[-1])
        lines[k][2].append(m[-1]-sq[-1])
        lines[k][3].append(m[-1]+sq[-1])
    
    names = ["../cache/braninnoise/PESVS_"+str(int(100*sp.log10(s)))+"_"+str(pwr)+"_"+str(i)+".p" for i in xrange(nreps)]
    results = search.multiPESVS(ojf,lb,ub,kernel,s,bd,lambda x,s:cfn(s),-7,-4,names)
    zr = [sp.log10(r[11].flatten()-braninmin) for r in results]
    
    m = sp.mean(sp.log10(Z),axis=0)
    C = [r[5] for r in results]
    for i in xrange(len(zr)):
        a[k+1].plot([sum(C[i][:j]) for j in xrange(len(C[i]))],zr[i],'b')
    pline[0].append(pwr)
    m_ = sp.mean([h[-1] for h in zr])
    pline[1].append(m_)
    sq = sp.sqrt(sp.var([h[-1] for h in zr]))
    pline[2].append(m_-sq)
    pline[3].append(m_+sq)
    
f_,a_ = plt.subplots(1)
for i in xrange(len(slist)):
    a_.plot(lines[i][0] ,lines[i][1])
    a_.fill_between(lines[i][0],lines[i][2],lines[i][3],facecolor='lightgreen',edgecolor='lightgreen',alpha=0.5)
a_.plot(pline[0],pline[1],'rx--')
a_.fill_between(pline[0],pline[2],pline[3],facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
f.savefig("../figs/costpower.png")
plt.show()