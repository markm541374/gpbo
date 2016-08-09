#!/usr/bin/env python2
#encoding: UTF-8

# Fitting hyperparameters for some power data, not connected to anything else yet

import scipy as sp
import numpy.random as npr
import pandas as pd
from matplotlib import pyplot as plt
#import seaborn as sns
import time
import gpbo.core.GPdc
import gpbo.core.search
import gpbo.core.OPTutils

days = 32
t0=time.clock()
df = pd.read_csv('../data/DemandData_Historic-2015.csv')
t1=time.clock()
print 'read time {0:e}'.format(t1-t0)
N = df.shape[0]
n = min(N,days*48)

dlb = 0.
dub = float(days)
print '{0:d} datapoints'.format(n)
X = sp.array([df.index.values[:n]]).T/48.
Y = sp.array([df.indo.values[:n]]).T/1000.
offs  = sp.mean(Y)
Y-=offs
f,a = plt.subplots(3)

f2,a2=plt.subplots(1)
a = sp.hstack([a,a2])

a[0].plot(X,Y,'g.')

S = sp.ones([n,1])*1e-6
D = [[sp.NaN]]*n

pm = sp.array([1.,-1.])
ps = sp.array([1.,1.])

def ojf(x,s,d,override=False):
    #print "called ojf: "+str(x)
    try:
        x=x.flatten(0)
    except:
        pass
    
    xlow = [-2.,-2.]
    xupp = [2.,2.]
    
    xthis = [xlow[i]+0.5*(xin+1)*(xupp[i]-xlow[i]) for i,xin in enumerate(x)]
    hyp = [10**i for i in xthis]
    
    print hyp
    t0=time.clock()
    llk = sp.clip(
        gpbo.core.GPdc.GP_LKonly(X, Y, S, D, gpbo.core.GPdc.kernel(gpbo.core.GPdc.MAT52, 1, sp.array(hyp))).plk(pm, ps), -1e60, 1e60)
    
    t1=time.clock()
    if llk<-1.:
        out = sp.log(-llk)+1.
    else:
        out = -llk
    print "--->llk: {0} {1}    t: {2}".format(llk,out,t1-t0)
    
    return [out,t1-t0]

def ojfa(x,s,d,override=False):
    #print "called ojf: "+str(x)
    
    try:
        x=x.flatten(0)
    except:
        pass
    
    xlow = [-2.,-2.]
    xupp = [2.,2.]
    
    xthis = [xlow[i]+0.5*(xin+1)*(xupp[i]-xlow[i]) for i,xin in enumerate(x[1:])]
    hyp = [10**i for i in xthis]
    #hyp = [10**i for i in x.flatten()[1:]]
    
    
    print hyp
    t0=time.clock()
    sub = x.flatten()[0]
    npts = int((1.-0.9*sub)*n)
    if override:
        npts=n
    print "subsampling {0} of {1} at x[0]={2}".format(npts,n,x.flatten()[0])
    ps = npr.choice(range(n),size=npts, replace=False)
    Xd = sp.vstack([X[i] for i in ps])
    Yd = sp.vstack([Y[i] for i in ps])
    Sd = sp.vstack([S[i] for i in ps])
    Dd = [[sp.NaN]]*npts
    
    llk = gpbo.core.GPdc.GP_LKonly(Xd, Yd, Sd, Dd, gpbo.core.GPdc.kernel(gpbo.core.GPdc.MAT52, 1, sp.array(hyp))).plk(pm, ps)
    t1=time.clock()
    if llk<-1.:
        out = sp.log(-llk)+1.
    else:
        out = -llk
    print "--->llk: {0} {1}    t: {2}".format(llk,out,t1-t0)
    
    return [out,t1-t0]


d=2
kindex = gpbo.core.GPdc.MAT52CS
prior = sp.array([0.]+[-1.]*d+[-2.])
sprior = sp.array([1.]+[1.]*d+[2.])
kernel = [kindex,prior,sprior]

#lets start with EI

lb = sp.array([[-1]*d])
ub = sp.array([[1]*d])

budget = 20
fnames = ['../cache/fith/EI{}.p'.format(i) for i in xrange(5)]
statesEI= gpbo.core.search.multiMLEFS(ojf, lb, ub, kernel, 1., budget, fnames)


fnames = ['../cache/fith/PE{}.p'.format(i) for i in xrange(5)]
statesPE= gpbo.core.search.multiPESFS(ojf, lb, ub, kernel, 1., budget, fnames)


kindex = gpbo.core.GPdc.MAT52CS
prior = sp.array([0.]+[-1.]*(d+1)+[-2.])
sprior = sp.array([1.]*(d+2)+[2.])
kernel = [kindex,prior,sprior]

fnames = ["../cache/fith/PI{}.p".format(i) for i in xrange(5)]
statesPI = gpbo.core.search.multiPESIPS(ojfa, lb, ub, kernel, 10, fnames)


x=[]
y=[]
for stateEI in statesEI:
    a[2].plot(stateEI[11],'r')
    #a[3].plot([sum(stateEI[5][:i]) for i in xrange(len(stateEI[5]))],stateEI[11],'r')
    x.append([sum(stateEI[5][:i]) for i in xrange(len(stateEI[5]))])
    y.append(stateEI[11].flatten())
    print 'reccomended under EI: {} : {}'.format([10**i for i in stateEI[4][-1]],ojf(stateEI[4][-1],None,None)[0])
X_,Y_,lb_,ub_ = gpbo.core.OPTutils.mergelines(x, y)
a[3].fill_between(X_,lb_,ub_,facecolor='lightcoral',edgecolor='lightcoral',alpha=0.5)
a[3].plot(X_,Y_,'r')

x=[]
y=[]
for statePE in statesPE:
    a[2].plot(statePE[11],'b')
    #a[3].plot([sum(statePE[5][:i]) for i in xrange(len(statePE[5]))],statePE[11],'b')
    x.append([sum(statePE[5][:i]) for i in xrange(len(statePE[5]))])
    y.append(statePE[11].flatten())
    print 'reccomended under PE: {} : {}'.format([10**i for i in statePE[4][-1]],ojf(statePE[4][-1],None,None)[0])
    
X_,Y_,lb_,ub_ = gpbo.core.OPTutils.mergelines(x, y)
a[3].fill_between(X_,lb_,ub_,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
a[3].plot(X_,Y_,'b')


x=[]
y=[]
for statePI in statesPI:
    a[2].plot(statePI[11],'c')
    #a[3].plot([sum(statePI[5][:i]) for i in xrange(len(statePI[5]))],statePI[11],'c')
    x.append([sum(statePI[5][:i]) for i in xrange(len(statePI[5]))])
    y.append(statePI[11].flatten())

    print 'reccomended under PI: {} : {}'.format([10**i for i in statePI[4][-1][1:]],ojf(statePI[4][-1][1:],None,None)[0])

X_,Y_,lb_,ub_ = gpbo.core.OPTutils.mergelines(x, y)
a[3].fill_between(X_,lb_,ub_,facecolor='lightgreen',edgecolor='lightgreen',alpha=0.5)
a[3].plot(X_,Y_,'g')
a[3].set_xscale('log')
a[3].set_xlabel('accumulated evaluation time(s)')
a[3].set_ylabel('GP negative log-likelihood')

    
#plt.plot(statePI[8],'g')
a[3].axis([0.1,10,8,24])
plt.show()