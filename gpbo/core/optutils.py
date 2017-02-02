# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
#import pyximport
#from numpy import get_include
#pyximport.install(setup_args={'include_dirs': get_include()})
from __future__ import print_function
xrange=range
import numpy as np
import scipy as sp
from scipy.optimize import minimize
from scipy.stats import gamma
import sys
import os
from scipy import linalg as spl
import time
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from gpbo.core import GPdc
from matplotlib import pyplot as plt
import DIRECT
from gpbo.core import ESutils
def cosines(x,s,d):
    x.resize([1,x.size])
    assert(d==[sp.NaN])
    
    u = 1.6*x[0,0]-0.5
    v = 1.6*x[0,1]-0.5
    f = 1.-(u**2 + v**2 -0.3*sp.cos(3*sp.pi*u)-0.3*sp.cos(3*sp.pi*v)+0.7)
    
    if s==0.:
        noise = 0.
    else:
        noise = sp.random.normal(scale=sp.sqrt(s))
    return [f +noise,1.]

def digitfn():
    from sklearn import datasets, svm
    digits = datasets.load_digits()
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))
    def digitojf(x,s,d,override=False):
        x.resize([1,x.size])
        g = 10**((x[0,0]-1.)*3)
        c = 6.*x[0,1]
        # Create a classifier: a support vector classifier
        classifier = svm.SVC(gamma=g, coef0=c,kernel='sigmoid')
        # We learn the digits on the first half of the digits
        classifier.fit(data[:n_samples / 2], digits.target[:n_samples / 2])
        # Now predict the value of the digit on the second half:
        r = classifier.score(data[n_samples / 2:],digits.target[n_samples / 2:])
        return [1-r,1.]
    return digitojf

def quad(x,s,d):
    assert(d==[sp.NaN])
    f = sum((x.flatten()-0.1)**2)
    if s==0.:
        noise = 0.
    else:
        noise = sp.random.normal(scale=sp.sqrt(s))
    return [f +noise,1.]

bananamin = sp.array([0.2,0.2])
def genbanana(ignores=-1.,cfn = lambda x:1.,):
    def banana(x,s,d, override=False):
        
        assert(d==[sp.NaN])
        x.resize([1,x.size])
        u = 5.*x[0,0]
        v = 5.*x[0,1]
        a=1.
        b=100.
        f = 1e-3*((a-u)**2 + b*(v-u**2)**2)
        if ignores>0:
            s=ignores
        if s==0.:
            noise = 0.
            
        else:
            
            noise = sp.random.normal(scale=sp.sqrt(s))
        try:
            return [f+noise,cfn(s)]
        except ZeroDivisionError:
            return [f+noise,sp.inf]
    return banana

def genbranin(ignores=-1.,cfn = lambda x:1.,):
    def branin(x,s,d, override=False):
        
        assert(d==[sp.NaN])
        x.resize([1,x.size])
        u = x[0,0]*7.5 + 2.5
        v = x[0,1]*7.5 + 2.5
        
        f = (-1.275*(u/sp.pi)**2+5*u/sp.pi+v-6)**2 +(10.-5./(4*sp.pi))*sp.cos(u) + 10.
        if ignores>0:
            s=ignores
        if s==0.:
            noise = 0.
            
        else:
            
            noise = sp.random.normal(scale=sp.sqrt(s))
        try:
            return [f+noise,cfn(s)]
        except ZeroDivisionError:
            return [f+noise,sp.inf]
    return branin

def gencamel(ignores=-1.,cfn = lambda x:1.,):
    def camel(x,s,d, override=False):
        
        assert(d==[sp.NaN])
        x.resize([1,x.size])
        u = x[0,0]*5
        v = x[0,1]*7.5
        
        f = 4*u**2+u*v-4*v**2-2.1*u**4+4*v**4+(u**6)/3.
        if ignores>0:
            s=ignores
        if s==0.:
            noise = 0.
            
        else:
            
            noise = sp.random.normal(scale=sp.sqrt(s))
        try:
            return [f+noise,cfn(s)]
        except ZeroDivisionError:
            return [f+noise,sp.inf]
    return camel


def gensquexpdraw(d,lb,ub,ignores=-1):
    nt=14
    [X,Y,S,D] = ESutils.gen_dataset(nt, d, lb, ub, GPdc.SQUEXP, sp.array([1.5] + [0.30] * d))
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, d, sp.array([1.5] + [0.30] * d)))
    def obj(x,s,d,override=False):
        #print [x,s,d]
        if ignores>0:
            s=ignores
        if s==0. or override:
            noise = 0.
        else:
            noise = sp.random.normal(scale=sp.sqrt(s))
        print( "EVAL WITH NOISE: "+str(noise) + "FROM S= "+str(s))
        return [G.infer_m(x,[d])[0,0]+noise,1.]
    def dirwrap(x,y):
        z = G.infer_m(x,[[sp.NaN]])[0,0]
        #z = obj(x,0.,[sp.NaN])
        return (z,0)
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    
    return [obj,xmin,ymin]

def gensquexpIPdraw(d,lb,ub,sl,su,sfn,sls,cfn):
    #axis = 0 value = sl
    #d dimensional objective +1 for s
    nt=25
    #print sp.hstack([sp.array([[sl]]),lb])
    #print sp.hstack([sp.array([[su]]),ub])
    [X,Y,S,D] = ESutils.gen_dataset(nt, d + 1, sp.hstack([sp.array([[sl]]),lb]).flatten(), sp.hstack([sp.array([[su]]),ub]).flatten(), GPdc.SQUEXP, sp.array([1.5] + [sls] + [0.30] * d))
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, d + 1, sp.array([1.5] + [sls] + [0.30] * d)))
    def obj(x,s,d,override=False):
        x = x.flatten()
        if sfn(x)==0. or override:
            noise = 0.
        else:
            noise = sp.random.normal(scale=sp.sqrt(sfn(x)))
        
        return [G.infer_m(x,[d])[0,0]+noise,cfn(x)]
    def dirwrap(x,y):
        z = obj(sp.array([[sl]+[i for i in x]]),sl,[sp.NaN],override=True)
        return (z,0)
    [xmin0,ymin0,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    lb2 = xmin0-sp.ones(d)*1e-4
    ub2 = xmin0+sp.ones(d)*1e-4
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb2,ub2,user_data=[], algmethod=1, maxf=89000, logfilename='/dev/null')
    #print "RRRRR"+str([xmin0,xmin,ymin0,ymin,xmin0-xmin,ymin0-ymin])
    return [obj,xmin,ymin]






def bounds(Xs,Ys,ns=100):
    #use a gp to infer mean and bounds on sets of x/y data that have diffent x
    #f,a = plt.subplots(2)
    #for i in xrange(len(Ys)):
    #    a[0].plot(Xs[i],Ys[i])
    
    X = sp.hstack(Xs)
    np = X.size
    Y = sp.hstack(Ys)
    X.resize([np,1])
    Y.resize([np,1])
    #a[1].plot(X,Y,'r.')
    np = X.size
    S = sp.zeros(np)
    D = [[sp.NaN]]*np
    ki = GPdc.MAT52CS
    mprior = sp.array([1.,2.,1.])
    sprior = sp.array([2.,2.,2.])
    #MAPH = GPdc.searchMAPhyp(X,Y,S,D,mprior,sprior, ki,mx=500)
    MAPH = sp.array([0.5,5.,0.3])
    g = GPdc.GPcore(X, Y, S, D, GPdc.kernel(ki, 1, MAPH))
    sup = sp.linspace(min(X),max(X),ns)
    [m,V] = g.infer_diag_post(sup,[[sp.NaN]]*ns)
    std = sp.sqrt(V+MAPH[2])
    #plt.fill_between(sup.flatten(),(m-std).flatten(),(m+std).flatten(),facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
    #a[1].plot(sup,m.flatten(),'b')
    return [sup,m,std]

from scipy.interpolate import interp1d
    
def mergelines(x,y):
    minx = max([min(i) for i in x])
    maxx = min([max(i) for i in x])
    fs = []
    for i in xrange(len(x)):
        #print [x[i].shape,y[i].shape]
        
        fs.append(interp1d(x[i],y[i]))
    X = [i for i in sorted(sp.hstack([sp.array(j) for j in x])) if i<=maxx and i>=minx]
    np = len(X)
    X=sp.array(X)
    Y=sp.empty(np)
    ub=sp.empty(np)
    lb=sp.empty(np)
    for i in xrange(np):
        q = [j(X[i]) for j in fs]
        Y[i] = sp.mean(q)
        v = sp.var(q)
        ub[i] = Y[i]+2.*sp.sqrt(v)
        lb[i] = Y[i]-2.*sp.sqrt(v)
    
    return X,Y,lb,ub

def silentdirect(f,l,u,*args,**kwargs):
    print( 'searching...')
    fileno = sys.stdout.fileno()
    with os.fdopen(os.dup(fileno), 'wb') as stdout:
        with os.fdopen(os.open(os.devnull, os.O_WRONLY), 'wb') as devnull:
            sys.stdout.flush();
            os.dup2(devnull.fileno(), fileno)  # redirect
            [xmin, ymin, ierror] = DIRECT.solve(f,l,u,*args,**kwargs)
        sys.stdout.flush();
        os.dup2(stdout.fileno(), fileno)
    print( 'direct found {} at {} {}'.format(ymin,xmin,ierror))
    return xmin,ymin,ierror

def overheadregression(x,y):
    x=list(x)
    y=list(y)
    N=len(x)

    def fun(x, theta):
        a = abs(theta[0])
        b = abs(theta[1])+1.
        c = abs(theta[2])
        return a * x ** b + c

    def llk(x, y, theta):
        x = list(x)
        y = list(y)
        N = len(x)
        t = abs(theta[3])
        acc = 0.
        for i in range(N):
            acc += -0.5 * ((fun(x[i], theta) - y[i]) ** 2)
        # print acc
        return acc / (t ** 2) - N * np.log(t)

    def lpr(theta, prshp, prscl):
        acc = 0.
        for i in range(4):
            # print(gamma.logpdf(theta[i],prshp[i],scale=prscl[i]),theta[i],prshp[i],prscl[i])
            acc += gamma.logpdf(abs(theta[i]), prshp[i], scale=prscl[i])
        # print acc
        return acc

    def lpr(theta, prmean, prstd):
        acc = 0.
        for i in range(4):
            acc += -0.5 * ((theta[i] - prmean[i]) ** 2) / prstd[i] ** 2
        return acc

    prshp = [0.25, 0.4, 5., 0.5]
    prscl = [0.4,0.05, 5., 1.]

    def f(theta):
        v = -llk(x, y, theta) - lpr(theta, prshp, prscl)
        # print theta,v
        return v

    res = minimize(f,[1.,1.5,1.,0.1])
    para = [abs(i) for i in res.x]
    #print( x,y)
    #from matplotlib import pyplot as plt
    z= sp.linspace(0,50,200)
    #plt.figure()
    #plt.plot(x,y,'r.')
    #plt.plot(z,map(lambda q:fun(q,para),z),'b')
    #plt.show()
    #plt.savefig('/home/mark/Desktop/zzz{}.png'.format(N))
    #plt.close()
    para[1]+=1
    print( 'overheadfit  {}*n^{} +{}+ N(0,{})'.format(*para))
    return para

def geteffectiveoverhead(optstate,nrandinit):
    coefs = overheadregression(range(1, optstate.n + 1)[nrandinit:], optstate.aqtime[nrandinit:])
    print('remaining {}'.format(optstate.remaining))
    print('CEVavg = {}'.format(sp.mean(optstate.c)))
    remain = optstate.remaining
    cev = sp.mean(optstate.c)
    over = lambda i: coefs[0] * i ** coefs[1] + coefs[2]
    lastover = optstate.aqtime[-1]

    # remaining steps under constant overhead
    R = int(remain / (lastover + cev))
    # remaining steps under growth prediction
    Nr = 0
    acc = 0.
    while remain > acc + Nr * cev:
        Nr += 1
        acc += over(Nr + optstate.n)
    print('const remain {} var remain {}'.format(R, Nr))
    print('under const acq: {} ev: {}'.format(R * lastover, R * cev))
    print('under predi acq: {} ev: {}'.format(acc, Nr * cev))
    print('current overhead: {} predicted av overhead: {}'.format(lastover, acc / float(Nr)))
    return acc/float(Nr)