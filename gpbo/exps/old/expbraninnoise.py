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
print 'start'
print 'start'
d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])
pwr = 0.1
def cfn(s):
    ##print s
    #print 'cfn'+str(((1e-6)/s)**pwr)
    return ((1e-6)/s)**pwr
ojf = gpbo.core.optutils.genbranin(cfn=cfn)
braninmin = 0.39788735772973816
kindex = gpbo.core.GPdc.MAT52
prior = sp.array([0.]+[-1.]*d)
sprior = sp.array([1.]*(d+1))
kernel = [kindex,prior,sprior]
nreps = 2
bd = 3.
slist = [1e3,1e-3]
print 'start'
f,a_ = plt.subplots(2)
fnull,anull=plt.subplots(1)
a = sp.hstack([anull,a_])
import os
tmp=[]
fc = ['lightblue','lightgreen']
lc = ['blue','green']
for ii,s in enumerate(slist):
    
    #names = ["../cache/braninnoise/EIMLE_"+str(int(100*sp.log10(s)))+"_"+str(pwr)+"_"+str(i)+".p" for i in xrange(nreps)]
    #results = search.multiMLEFS(ojf,lb,ub,kernel,s,bd,names)
    #yr = [r[11].flatten() for r in results]
    #C = results[0][5]
    #print 'start'
    b_ = min(bd,200.*cfn(s))
    names = ["../cache/braninnoise/PESFS_"+str(int(100*sp.log10(s)))+"_"+str(pwr)+"_"+str(i)+".p" for i in xrange(nreps)]
    results = gpbo.core.search.multiPESFS(ojf, lb, ub, kernel, s, b_, names)
    zr = [r[11].flatten() for r in results]
    C = results[0][5]
    tmp.append([ len(r[5]) for r in results])
        
    #Z = sp.vstack(yr)-braninmin
    #m = sp.mean(sp.log10(Z),axis=0)
    #v = sp.var(sp.log10(Z),axis=0)
    
    #sq = sp.sqrt(v)
    #print len(C)
    #print m.size
    #a[2].fill_between(sp.array([sum(C[:j]) for j in xrange(len(C))]),(m-sq).flatten(),(m+sq).flatten(),facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
    #a[2].plot([sum(C[:j]) for j in xrange(len(C))],m.flatten(),'x-')
    
    
    Z = sp.vstack(zr)-braninmin
    m = sp.mean(sp.log10(Z),axis=0)
    v = sp.var(sp.log10(Z),axis=0)
    sq = sp.sqrt(v)
    a[2].fill_between(sp.array([sum(C[:j]) for j in xrange(len(C))]),(m-sq).flatten(),(m+sq).flatten(),facecolor=fc[ii],edgecolor=fc[ii],alpha=0.5)
    a[2].plot([sum(C[:j]) for j in xrange(len(C))],m.flatten(),'.-',color=lc[ii])
    
    
    
    print 'start'
    f.savefig("tmp.png")
    

s=1e-1

#for i in xrange(nreps):
#    fname = "../cache/rosennoise/PESVS_"+str(int(100*sp.log10(s)))+"_"+str(pwr)+"_"+str(i)+".p"
#    [X,Y,S,D,R,C2,T,Tr,Ymin,Xmin,Yreg, Rreg] = search.PESVS(ojf,lb,ub,kernel,s,ba,lambda x,s:cfn(s),-9,-1,fname)
lsu = 3
lsl = -7
names = ["../cache/braninnoise/PESVS_"+str(int(100*sp.log10(s)))+"_"+str(pwr)+"_"+str(i)+"_"+str(lsl)+"_"+str(lsu)+".p" for i in xrange(nreps)]
results = gpbo.core.search.multiPESVS(ojf, lb, ub, kernel, s, bd, lambda x, s:cfn(s), lsl, lsu, names)
Rz = [sp.log10(sp.array(r[11])-braninmin).flatten() for r in results]
Cz = [[sum(r[5][:j]) for j in xrange(len(r[5]))] for r in results]
#for i in xrange(nreps):
#    a[2].plot(Cz[i],Rz[i].flatten(),'rx-')
#[a[1].semilogy([1e-6/(c**(1./pwr)) for c in r[5]],'cyan') for r in results]
x=[]
y=[]
for r in results:
    y.append(sp.array(sp.log10([1e-6/(c**(1./pwr)) for c in r[5]])))
    print y
    x.append(range(len(r[5])))
X_,Y_,lb_,ub_ = gpbo.core.optutils.mergelines(x, y)
a[1].fill_between(X_,lb_,ub_,facecolor='lightcoral',edgecolor='lightcoral',alpha=0.5)
a[1].plot(X_,Y_,'r')
    #a[2].plot(sp.array([sum(C2[:j]) for j in xrange(len(C2))]).flatten(),(sp.log(Yreg)).flatten(),'ro-')


def linint(i,c):
    if c>Cz[i][-1]:
        return sp.NaN
    else:
        j=0
        while Cz[i][j]<c:
            j+=1
    return Rz[i][j-1]+(c-Cz[i][j-1])*(Rz[i][j]-Rz[i][j-1])/(Cz[i][j]-Cz[i][j-1])

#for i in xrange(nreps):
#    a[2].plot(Cz[i],Rz[i],'c',alpha=0.5)
    
li = sp.empty(1000)
vi = sp.empty(1000)
for i,c in enumerate(sp.linspace(0,bd,1000)):
    tmp=[]
    for j in xrange(nreps):
        
        t2 = linint(j,c)
        if not sp.isnan(t2):
            tmp.append(t2)
    if len(tmp)<=1:
        li[i]=sp.NaN
        vi[i]=sp.NaN
    else:
        li[i]= sp.mean(tmp)
        vi[i]= sp.sqrt(sp.var(tmp))
a[2].fill_between(sp.linspace(0,bd,1000),li-vi,li+vi,facecolor='salmon',edgecolor='salmon',alpha=0.5)
a[2].plot(sp.linspace(0,bd,1000),li,'r')

#[sup,m,sd]=OPTutils.bounds(Cz,Rz)
#a[2].fill_between(sup.flatten(),(m-sd).flatten(),(m+sd).flatten(),facecolor='salmon',edgecolor='salmon',alpha=0.5)
#a[2].plot(sup,m.flatten(),'darkred')


print [len(r[5]) for r in results]
print tmp
sx = sp.logspace(0,-8,100)
cx = map(cfn,sx)
a[0].loglog(sx,cx)
if pwr==0.4:
    a[2].axis([0,0.5,-4,2])
a[2].set_xlabel('accumulated evaluation cost')
a[1].set_xlabel('steps')
a[2].set_ylabel('log10 regret')
a[1].set_ylabel('log10 observation variance')
f.savefig("../figs/braninnoise.png")
plt.show()