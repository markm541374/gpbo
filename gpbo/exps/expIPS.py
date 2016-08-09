#!/usr/bin/env python2
#encoding: UTF-8

#Run the env variable PES on synthetic functions, compare regular PES fixed at s=0, 0.25, 1. The !=0 should converge to a high regret according to the min in the plane they are in, the =0 and env var should kee[p improving


import gpbo.core.OPTutils
import scipy as sp
from matplotlib import pyplot as plt
import gpbo.core.GPdc
import gpbo.core.ESutils
import gpbo.core.search
from tqdm import tqdm, tqdm_gui
import DIRECT


seed=123
sp.random.seed(seed)
d=2
lb = sp.array([[-1.]*d])
ub = sp.array([[1.]*d])
sl=0.
su=1.
sfn = lambda x:1e-8
fls = 1.0
sls = (su-sl)*fls
dcc=1.0
cfn = lambda x: sp.exp(-dcc*x.flatten()[0])
[ojf,truexmin,ymin] = gpbo.core.OPTutils.gensquexpIPdraw(d, lb, ub, sl, su, sfn, sls, cfn)

#what are the best mins in planes away from true?
def planemin(xp):
    def dirwrap(x,y):
        z,c = ojf(sp.hstack([xp,x]),-1,[sp.NaN],override=True)
        return (z,0)
    [xm,ym,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=80000, logfilename='/dev/null')
    ye,c = ojf(sp.hstack([0.,xm]),-1,[sp.NaN],override=True)
    r=ye-ymin
    return [xm,ym,ye,r]

plane00 = planemin(0.)
plane025 = planemin(0.25)
plane05 = planemin(0.5)
plane001 = planemin(0.01)
plane10 = planemin(1.)
nreps=6
bd=35

kindex = gpbo.core.GPdc.MAT52CS
prior = sp.array([0.]+[-1.]*(d+1)+[-2.])
sprior = sp.array([1.]*(d+2)+[2.])
kernel = [kindex,prior,sprior]

names = ["../cache/IPS_/PESIPS_"+str(dcc)+"_"+str(fls)+"_"+str(seed)+"_"+str(i)+".p" for i in xrange(nreps)]
results = gpbo.core.search.multiPESIPS(ojf, lb, ub, kernel, bd, names)

f,a = plt.subplots(2)
aot = a[0].twinx()
Xa = [sp.array(r[0][:,0]) for r in results]
#Ca = [[sum(r[5][:j]) for j in xrange(len(r[5]))] for r in results]

x=[]
y=[]
w=[]

for r in results:
    #a[0].plot(r[0][:,0].flatten(),'b')
    
    #a[0].plot(r[5],'r')
    #a[1].plot((r[11].flatten()-ymin),'b')
    y.append(sp.array(sp.log(r[11].flatten()-ymin)))
    x.append(range(len(y[-1])))
    #a[0].plot([sum(r[5][:j]) for j in xrange(len(r[5]))],(r[11].flatten()-ymin),'b')
    w.append([sum(r[5][:j]) for j in xrange(len(r[5]))])

X_,Y_,lb_,ub_ = gpbo.core.OPTutils.mergelines(w, y)
a[0].fill_between(X_,lb_,ub_,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
a[0].plot(X_,Y_,'b')

X_,Y_,lb_,ub_ = gpbo.core.OPTutils.mergelines(x, y)
a[1].fill_between(X_,lb_,ub_,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
a[1].plot(X_,Y_,'b')

a[1].set_ylabel("log10 regret")
a[0].set_ylabel("log10 regret")
a[0].set_xlabel("accumulated cost")
a[1].set_xlabel("steps")
#a[1].set_xscale("log")
#a[0].set_xscale("log")
subset = [0.,0.25,1.]
c = ['g','r','c']

print "XXX"
print [truexmin,ymin]
print plane00
#a[2].plot([0,bd],[sp.log10(plane00[3])]*2,c[0],linestyle='--')
print plane025
#a[2].plot([1.,bd],[sp.log10(plane025[3])]*2,color=c[1],linestyle='--')
print plane10
#a[2].plot([1.,bd],[sp.log10(plane10[3])]*2,color=c[2],linestyle='--')
print plane001
#a[2].plot([1.,bd],[sp.log10(plane001[3])]*2,color='k',linestyle='--')
print "XXX"

kindex = gpbo.core.GPdc.MAT52CS
priora = sp.array([0.]+[-1.]*(d)+[-2.])
spriora = sp.array([1.]*(d+1)+[2.])
kernela = [kindex,priora,spriora]

for i,xs in enumerate(subset):
    def ojfa(x,s,d,override=False):
        return ojf(sp.hstack([[xs],x.flatten()]),s,d,override=override)
    names = ["../cache/IPS_/PESIS_"+str(xs)+"_"+str(dcc)+"_"+str(fls)+"_"+str(seed)+"_"+str(k)+".p" for k in xrange(nreps)]
    results = gpbo.core.search.multiPESIS(ojfa, lb, ub, kernela, bd, names)
    Rg = sp.log10(sp.vstack([(sp.array([ojf(sp.hstack([0.,r[4][j,:]]) ,0.,[sp.NaN],override=True)[0] for j in xrange(r[4].shape[0])])-ymin) for r in results]))
    mr = sp.mean(Rg,axis=0)
    sr = sp.sqrt(sp.var(Rg,axis=0))
    cacc = [sum(r[5][:j]) for j in xrange(len(r[5]))]
    a[0].fill_between(cacc,(mr-sr).flatten(),(mr+sr).flatten(),facecolor=c[i],edgecolor=c[i],alpha=0.5)
    a[0].plot(cacc,mr.flatten(),color=c[i])
    x=[]
    y=[]
    for r in results:
        #a[0].plot(r[5],color=c[i])
        z=sp.array([ojf(sp.hstack([0.,r[4][j,:]]) ,0.,[sp.NaN],override=True)[0] for j in xrange(r[4].shape[0])])
        #a[1].plot(sp.log10((z-ymin)),color=c[i])
        y.append(sp.log10((z-ymin)))
        x.append(range(len(y[-1])))
    X_,Y_,lb_,ub_ = gpbo.core.OPTutils.mergelines(x, y)
    a[1].fill_between(X_,lb_,ub_,facecolor=c[i],edgecolor=c[i],alpha=0.5)
    a[1].plot(X_,Y_,c[i])
        
        #a[2].plot([sum(r[5][:j]) for j in xrange(len(r[5]))],sp.log10(z-ymin),color=c[i])
a[1].set_xscale('log')
a[0].set_xscale('log')  

plt.savefig('../figs/expIPS.png')
plt.show()
