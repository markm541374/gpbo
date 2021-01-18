# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

#cython: profile=True
from __future__ import print_function
xrange=range
import numpy as np
import scipy as sp
from scipy.optimize import minimize
from scipy.stats import gamma
from scipy import stats as sps
from scipy import integrate as spi
import sys
import os
from scipy import linalg as spl
import time
import tqdm
from gpbo.core import GPdc
from matplotlib import pyplot as plt
import DIRECT
from gpbo.core import ESutils
from libc.math cimport exp, sqrt, M_SQRT2, M_PI,erf,log


import logging
logger = logging.getLogger(__name__)
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
    t0=time.clock()
    try:
        fileno = sys.stdout.fileno()
    except:
        fileno = sys.stdout
    with os.fdopen(os.dup(fileno), 'wb') as stdout:
        with os.fdopen(os.open(os.devnull, os.O_WRONLY), 'wb') as devnull:
        #with os.fdopen(os.dup(fileno), 'wb') as devnull:
            sys.stdout.flush();
            os.dup2(devnull.fileno(), fileno)  # redirect
            [xmin, ymin, ierror] = DIRECT.solve(f,l,u,*args,**kwargs)
        sys.stdout.flush();
        os.dup2(stdout.fileno(), fileno)
    #print( 'direct found {} at {} {}'.format(ymin,xmin,ierror))
    print('directtime {}'.format(time.clock()-t0))
    return xmin,ymin,ierror


def noisydirect(f,l,u,*args,**kwargs):
    print( 'searching...')
    t0=time.clock()
    [xmin, ymin, ierror] = DIRECT.solve(f,l,u,*args,**kwargs)
    print('directtime {}'.format(time.clock()-t0))
    return xmin,ymin,ierror

def silentdirectwrapped(f,l,u,*args,**kwargs):
    def dwrap(x,aux):
        return f(x),0
    return silentdirect(dwrap,l,u,*args,**kwargs)

def noisydirectwrapped(f,l,u,*args,**kwargs):
    def dwrap(x,aux):
        return f(x),0
    return noisydirect(dwrap,l,u,*args,**kwargs)

def boundedlocal(f,l,u,x0,*args,**kwargs):
    d = len(u)
    res = minimize( f,x0,method='L-BFGS-B',bounds=tuple([(l[j],u[j]) for j in range(d)]),options=kwargs)
    xmin,ymin,ierror = res.x,res.fun,res.message
    return xmin,ymin,ierror

def twopartopt(f,l,u,dargs,largs):
    dxmin,dymin,dierror = silentdirectwrapped(f,l,u,**dargs)
    xmin,ymin,ierror = boundedlocal(f,l,u,dxmin,**largs)
    logger.debug('optresult {} at {}'.format(str(ymin), str(xmin)))
    return xmin,ymin,ierror






def multilocal(f,l,u,*args,**kwargs):
    print('searching multistart nmax={}'.format(kwargs['maxf']))
    D = len(l)
    start = sp.empty(D)
    fevacc = 100
    bounds=tuple([(l[j],u[j]) for j in range(D)])
    minf = sp.Inf
    Xinit = sp.empty([100,D])
    Finit = sp.empty(100)
    for i in xrange(100):
        for j in xrange(D):
            Xinit[i,j] = sp.random.uniform(l[j],u[j])
        Finit[i] = f(Xinit[i,:])

    while fevacc<kwargs['maxf']:
        k = sp.argmin(Finit)
        start = Xinit[k,:]
        Finit[k]=sp.Inf
        res = minimize(f ,start,method='L-BFGS-B',bounds=bounds,options={'ftol':0.0001,'maxfun':60+20*D})

        #print('from {} found {} {} {} {}'.format(start,res.x,res.fun,res.message,res.nfev))
        if res.fun<minf:
            minf=res.fun
            result = res.x
            message = res.message
        fevacc+=res.nfev
    return result,minf,message
def llk(X,Y,theta,pr_alpha, pr_beta):
    X=list(X)
    Y=list(Y)
    N=len(X)
    a,b,c,s = theta
    accf=-N*np.log(s)
    accD = np.array([0.,0.,0.,-N/s])
    for i in range(4):
        accf+= (pr_alpha[i]-1.)*np.log(theta[i])-theta[i]*pr_beta[i]
        accD[i] += (pr_alpha[i]-1.)/theta[i] - pr_beta[i]
    for i in xrange(N):
        error = (a*X[i]**b+c-Y[i])/(s**2)
        accD[0] -= error*(X[i]**b)
        accD[1] -= error*(a*np.log(X[i])*X[i]**b)
        accD[2] -= error
        accD[3] += (error**2)*s
        accf -= (0.5*error**2)*s**2
    #print theta,accf
    return accf,accD


def overheadregression(X,Y):
    #print(X,Y)
    pr_alpha = [4.,8.,10.,2.]
    pr_beta = [2.,4.,1,0.5]
    def f(x):
        l,g =llk(X,Y,x,pr_alpha,pr_beta)
        #print x,-l
        return -l,-g

    res = minimize(f,[1.,2.,1.,1.],method='L-BFGS-B',jac=True,bounds=((1e-6,None),(1.,None),(1e-6,None),(1e-6,None)))#,bounds=((0.,1e3),(0.,1e3),(0.,1e3),(0.,1e3)))
    para = list(map(abs,res.x))
    print( 'fit model {}*x^{} + {} +N(0,{})'.format(*para))
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

def gpplot(meanaxis,varaxis,G,lb,ub,ns=50,nc=20):
    x_ = sp.linspace(lb[0],ub[0], ns)
    y_ = sp.linspace(lb[1],ub[1], ns)
    z_ = sp.empty([ns, ns])
    s_ = sp.empty([ns, ns])
    for i in range(ns):
        for j in range(ns):
            m_, v_ = G.infer_diag_post(sp.array([y_[j], x_[i]]), [[sp.NaN]])
            z_[i, j] = m_[0, 0]
            s_[i, j] = sp.sqrt(v_[0, 0])
    CS = meanaxis.contour(x_, y_, z_, nc)
    meanaxis.clabel(CS, inline=1, fontsize=10)
    CS = varaxis.contour(x_, y_, s_, nc)
    varaxis.clabel(CS, inline=1, fontsize=10)
    return

def Hvec2H(Hvec,d):
    Hvec=Hvec.flatten()
    H = sp.empty(shape=[d,d])
    k = 0
    for i in xrange(d):
        for j in xrange(i + 1):
            H[i, j] = H[j, i] = Hvec[k]
            k += 1
    return H

def gpGH(G,x):
    """
    get the mean and joint covariance of gradient and hessian of a GP returns vectorized H after G
    :param G: a GP
    :param x: location X
    :return:
    """
    d = G.D

    divs = [[]] * ((d * (d + 1) / 2)+d)
    k = d
    for i in xrange(d):
        divs[i]=[i]
        for j in xrange(i + 1):
            divs[k] = [i, j]
            k += 1
    X = sp.vstack([x] * ((d*(d + 1)/2)+d) )

    M,varM = G.infer_full_post(X,divs)
    G = M[:,:d]
    varG = varM[:d,:d]
    Hvec = M[:,d:]
    varHvec = varM[d:,d:]
    H = Hvec2H(Hvec,d)
    return G,varG,H,Hvec,varHvec,M,varM
def allvec2all(L,d):

    y = L[0]
    M = L.flatten()[1:]
    G = M[:d]
    Hvec = M[d:]
    H = Hvec2H(Hvec,d)
    return y,G,H

def gpYGH(G,x):
    """
    get the mean and joint covariance of y gradient and hessian of a GP returns vectorized H after G after y
    :param G: a GP
    :param x: location X
    :return:
    """
    d = G.D

    divs = [[]] * ((d * (d + 1) / 2)+d+1)
    k = d
    divs[0]=[np.NaN]
    for i in xrange(d):
        divs[i+1]=[i]
        for j in xrange(i + 1):
            divs[k+1] = [i, j]
            k += 1
    X = sp.vstack([x] * ((d*(d + 1)/2)+d+1) )
    L,varL = G.infer_full_post(X,divs)
    y = L[0,0]
    vary = varL[0,0]
    M = L[:,1:]
    varM = L[1:,1:]
    G = M[:,:d]
    varG = varM[:d,:d]
    Hvec = M[:,d:]
    varHvec = varM[d:,d:]
    H = Hvec2H(Hvec,d)
    return y,vary,G,varG,H,Hvec,varHvec,L,varL

def drawconditionH(G,varG,H,Hvec,varHvec,M,varM):
    d=G.size
    #Hdist = sp.stats.multivariate_normal(Hvec.flatten(), varHvec)
    H = GPdc.draw(Hvec.flatten(),varHvec,1).flatten()
    Khg = varM[:d,d:]
    #H = Hdist.rvs()
    choH = sp.linalg.cho_factor(varHvec)
    Gm = G + Khg.dot(spl.cho_solve(choH,H))
    Gv = varG - Khg.dot(spl.cho_solve(choH,Khg.T))
    Hd = Hvec2H(H,d)
    return Gm,Gv,Hd

def plotprobstatellipse(cG,H,x,ax,logr=False):
    # svd on cov of grad
    Ucg, Ecg, Vcg = spl.svd(cG)
    # new rotated covariance
    C0 = spl.solve(H, Ucg)
    varP = C0.dot(sp.diag(Ecg).dot(C0.T))
    # svd of new cov
    Uvp, Evp, Vvp = spl.svd(varP)
    circ = sp.empty([2, 200])
    for i in xrange(200):
        theta = 2. * sp.pi * i / 199.-sp.pi
        circ[:, i] = Uvp.dot(sp.array([sp.sin(theta) * 3*sp.sqrt(Evp[0]), sp.cos(theta) * sp.sqrt(Evp[1])])) + sp.array(
            [[j for j in x]])
    if logr:
        plotaslogrtheta(circ[0,:],circ[1,:],x[0],x[1],ax,'g')
    else:
        ax.plot(circ[0, :], circ[1, :], 'g')
    return


def probgppve(G,x,tol=1e-3,dropdims=[]):
    nsam = 10*int(1./tol)+1
    Gr, varG, H, Hvec, varHvec, M, varM = gpGH(G,x)
    d=G.D
    vHdraws = GPdc.draw(Hvec.flatten(),varHvec,nsam)
    pvecount = 0
    for i in xrange(nsam):
        Hdraw = Hvec2H(vHdraws[i,:], d)
        for i in reversed(sorted(dropdims)):
            Hdraw = np.delete(Hdraw,i,axis=0)
            Hdraw = np.delete(Hdraw,i,axis=1)
        try:
            sp.linalg.cholesky(Hdraw)
            pvecount += 1
        except:
            pass
    return (pvecount+1)/ float(nsam+2)

def plotaslogrtheta(X,Y,x0,x1,ax,*args,**kwargs):
    R = sp.log10(sp.sqrt((X-x0)**2+(Y-x1)**2))
    T = sp.mod(sp.arctan((Y-x1)/(X-x0))+0.5*(1-sp.sign(X-x0))*sp.pi,2*sp.pi)
    p = T.argsort()
    return ax.plot(R[p],T[p],*args,**kwargs)

def drawpartitionmin(G,S,xm,rm,n):
    #distance to xmin
    xm = sp.array(xm)
    ns,d = S.shape
    #R is distance from xm for each S
    R = sp.empty(ns)
    for i in xrange(ns):
        R[i] = sp.linalg.norm(S[i,:]-xm)
    #O is indicies by distance from xm
    O = sp.argsort(R)
    split = sp.searchsorted(R[O],rm)+1
    S_ = sp.vstack([xm,S[O,:]])
    Z = G.draw_post(S_, [[sp.NaN]] *(ns+1),  n)
    Res = sp.empty([n,5])
    Res[:,1] = Z[:,:split].min(axis=1)
    Res[:,2] = Z[:,split:].min(axis=1)
    Res[:,3] = Z[:,0]
    Res[:,0] = Res[:,1:3].min(axis=1)
    Res[:,4] = Res[:,1:3].argmin(axis=1)

    argminin = Z[:,:split].argmin(axis=1)
    argminmax = argminin.max()

    maxRin = R[O[argminmax-1]]
    #print(str(argminin)+'\n'+str(argminmax)+'\n'+str(maxRin)+' ' +str(rm)+'\n'+str(R[O[argminin]]))
    return Res, maxRin


def rline(unitvec,rmax,condition,htarget=1e-6):
    """
    :param unitvec: direction from origin for the linesearch
    :param rmax: outer limit
    :param condition: boolean fn of x
    :param htarget: resolution for binary search
    :return: maximum radius where condition is True, numevals made
    """
    if condition(unitvec*rmax):
        return rmax,1
    left=0
    right=rmax
    n=0
    while (right-left)>htarget:
        r = 0.5*(left+right)
        c = int(condition(unitvec*r))
        n+=1
        left = (1-c)*left + c*r
        right = (1-c)*r + c*right
    return left,n

def ballradsearch(d,rmax,condition,ndirs=200,lineSh=1e-6):
    evcount=0
    R=rmax
    for i in tqdm.tqdm(range(ndirs)):
        x = sp.stats.norm.rvs(sp.zeros(d)) #draw an unnormalized vec
        xunit = x/sp.linalg.norm(x) # normalize
        r,n = rline(xunit,R,condition,htarget=lineSh) # linesearch up to current max
        R = min(R,r) #drop max to new result
        if R==0.:
            break
        evcount+=n
    return R

cdef double NP(double x, double mu=0., double sigma=1.):
    return (1./(sqrt(2*M_PI)*sigma))*exp(-0.5*((x-mu)/sigma)**2)

cdef double NC(double x, double mu=0., double sigma=1.):
    return 0.5*(1+erf((x-mu)/(sigma*M_SQRT2)))

cdef class bigaussmin(object):
    """
    min of two gaussians as described by
    https://www.gwern.net/docs/conscientiousness/2008-nadarajah.pdf
    """
    cdef double x1, x2,std1,std2,rho,alpha,integraltol
    cdef object doubledist

    def __init__(self,x1=1.,std1=1.,x2=1.,std2=1.,rho=0.):
        self.x1 = x1
        self.x2 = x2
        self.std1 = std1
        self.std2 = std2
        self.rho = rho
        cov=self.std1*self.std2*rho
        self.alpha  = 1./sqrt(1-self.rho**2)
        self.doubledist = spl.cholesky([[std1**2,cov],[cov,std2**2]]).T
        self.integraltol=1e-9

        return
    def __repr__(self):
        return 'bigaussian minimum with means {} {} vars {} {} and corr {}'.format(self.x1,self.x2,self.std1,self.std2,self.rho)

    cdef double f1(self,double y):
        cdef double t0 = NP((y-self.x1)/self.std1)/self.std1
        cdef double a = self.rho*(y-self.x1)*self.alpha/self.std1
        cdef double b = (y-self.x2)*self.alpha/self.std2
        cdef double t1 = NC(a-b)
        return t0*t1
    cdef double f2(self,double y):
        cdef double t0 = NP((y-self.x2)/self.std2)/self.std2
        cdef double a = self.rho*(y-self.x2)*self.alpha/self.std2
        cdef double b = (y-self.x1)*self.alpha/self.std1
        cdef double t1 = NC(a-b)
        return t0*t1

    cpdef double pdf(self,double x):
        return self.f1(x)+self.f2(x)

    cpdef double cdf(self,double x):
        I,err = spi.quad(self.pdf,-sp.inf,x,epsabs=self.integraltol)
        return I

    def rvs(self,n):
        return (self.doubledist.dot(sp.random.normal(size=[2,n]))+sp.array([[self.x1],[self.x2]])).min(axis=0)

    def fit(self,X,startx=[0.1,1.,0.1,1.1,0.]):
        startx[1]=sqrt(startx[1])
        startx[3]=sqrt(startx[3])
        opt2para = lambda h:[h[0],h[1]**2,h[2],h[3]**2,(1.-1e-6)*2*sp.arctan(h[4])/sp.pi]
        def lk(h):
            dist = bigaussmin(*opt2para(h))
            cdef double l=0.
            cdef int i
            for i in range(X.size):
                l+=log(dist.pdf(X[i]))
            return -l
        R = sp.optimize.minimize(lk,startx,method='BFGS',options={'maxiter':500})
        self.__init__(*opt2para(R.x))
        print(self)
        print(R.x)
        return self

    def ERleft(self,x):
        I, err = spi.quad(lambda y:self.pdf(y)*(x-y),-sp.inf,x,epsabs=self.integraltol)
        return I
