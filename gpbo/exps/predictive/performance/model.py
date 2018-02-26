import numpy as np
import scipy as sp
from scipy.special import gamma
from scipy.optimize import minimize
from matplotlib import pyplot as plt
import os
import gpbo
import tqdm
import pickle

def decay2model(X,x0,y0,m0,x1,y1,m1):
    a0 = (10**y0)/((10**x0)**m0)
    a1 = (10**y1)/((10**x1)**m1)
    #return np.log(a0*X**m0 + a1*X**m1+1)
    return a0*X**m0 + a1*X**m1

def flsmodel(X,theta):

    yinit = theta[0]
    grad0 = theta[1]
    grad1 = theta[4]
    c0 = theta[2]+X[:,1]*theta[3]
    #softmax=np.exp(theta[5])
    return decay2model(X[:,0],1,yinit,grad0,0,c0,grad1)


normlogpdf = lambda x,loc,scale: -0.5*((x-loc)/scale)**2-0.5*np.log(2*np.pi*scale**2)
gammalogpdf = lambda x,shape,scale: (shape-1)*np.log(x)-x/scale-np.log(gamma(shape))-shape*np.log(scale)#-shape*np.log(scale)-np.log(gamma(shape))
tlogpdf = lambda x,v,scale: -0.5*(v+1)*np.log(1+((x/scale)**2)/v) + np.log(gamma(0.5*(v+1))/(np.sqrt(v*np.pi)*gamma(0.5*v))) - np.log(scale)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

path = 'data/results24'
files = os.listdir(path)

speedsub=7
def readdata(ls,ax):
    cachefile = os.path.join(path,'cache','{}_{}.p'.format(ls,speedsub))
    if os.path.isfile(cachefile):
        print('using cache {}'.format(cachefile))
        X,Y,plotdata = pickle.load(open(cachefile,'r'))
        for i,noise in enumerate([-2,-3,-4,-5,-6,-7,-8]):
            P = plotdata[i]
            ax.plot(P[0], P[2], color=colors[i], linestyle='-', label='$\sigma^2={}$'.format(noise))
            ax.fill_between(P[0],P[3],P[1],edgecolor=colors[i], linestyle='-',facecolor=colors[i],lw=0.0,alpha=0.1)
        return X,Y
    print('readdata {}'.format(ls))
    base = 'eihyp_3_{}_'.format(ls)
    x=[]
    y=[]
    plotdata=[]
    for i,noise in tqdm.tqdm(enumerate([-2,-3,-4,-5,-6,-7,-8])):
        names = [f for f in files if f.startswith(base+str(1000*noise))]

        D = [gpbo.optimize.readoptdata(os.path.join(path,n)) for n in names]

        xi = np.hstack([D[k]['index'][10:] for k in range(len(D))])
        x.append(np.vstack([xi,noise*np.ones_like(xi)]).T[::speedsub,:])
        y.append(np.hstack([D[k]['trueyatxrecc'][10:] for k in range(len(D))]).reshape([-1,1])[::speedsub])
        #a0[0].plot(x[-1],y[-1],color=colors[i],marker='.',linestyle='none')
        P = gpbo.opts.getmvint([D[k]['index'] for k in range(len(D))],[D[k]['trueyatxrecc'] for k in range(len(D))],logy=True,nstd=1)
        plotdata.append(P)
        ax.plot(P[0], P[2], color=colors[i], linestyle='-', label='$\sigma^2={}$'.format(noise))
        ax.fill_between(P[0],P[3],P[1],edgecolor=colors[i], linestyle='-',facecolor=colors[i],lw=0.0,alpha=0.1)
    X = np.vstack(x)
    Y = np.vstack(y)
    pickle.dump([X,Y,plotdata],open(cachefile,'w'))
    return X,Y
#for ls in [100,200,300,500,700,900,1200]:

def perfmodel(X,ls,t):
    m,c = np.split(t,2)
    tl = m*np.log(ls) + c
    mu = flsmodel(X,tl)
    stdlog = tl[-1]
    return mu, stdlog

def llk(X,Y,L,T):
    Yp,std = perfmodel(X,L,T)
    Ypl = np.log10(Yp).reshape(-1,1)
    #Yh = flsmodel(X,theta[:5]).reshape(-1,1)
    e = Ypl - np.log10(Y)
    #L = -np.sum(normlogpdf(e,0.,theta[5]))
    #L= -np.sum(normlogpdf(Y,Yh,Yh*theta[5]))
    #R = gammalogpdf(Y,theta[5],Yh/theta[5])
    R =  normlogpdf(Ypl-np.log10(Y),0.,std)
    L = -np.nansum(R[np.isfinite(R)])
    if np.isinf(L):
        pass
    return L

#def perlsllk(X,Y):
#    print('fitmodel')
#    theta0 = [-0.,-6,-0.25,0.5,-0.5,1.]
#    bds = ((-np.Inf,np.Inf),(-np.Inf,-1e-9),(-np.Inf,np.Inf),(1e-9,np.Inf),
#           (-np.Inf,-1e-9),(1e-9,np.Inf))
#    f = lambda t:llk(X,Y,t)
#    res = minimize(f,theta0,method='L-BFGS-B',bounds=bds)
#    print('{}'.format(res.x))
#    return res.x

L = [100,200,300,500,700,900,1200,1500]
def optfull():
    X=[]
    Y=[]
    print('readdata')
    F=[]
    A=[]
    for ls in L:
        f,a = plt.subplots()
        F.append(f)
        A.append(a)
        x,y = readdata(ls,a)
        X.append(x)
        Y.append(y)

    def fullllk(t):
        r=0
        #m,c = np.split(t,2)
        for i,ls in enumerate(L):
            #tl = m*np.log(ls/1000.) + c
            r+=llk(X[i],Y[i],ls/1000.,t)
        print(np.array_str(t,precision=3),r)
        return r

    print('fitmodel')
    theta0 = [0.,-1.,0.,0.,0.,0.,-0.,-6,-0.25,0.5,-0.5,1.]
    bds =  ((-np.Inf,np.Inf),)*len(theta0)
    print bds
    res = minimize(fullllk,theta0,method='L-BFGS-B',bounds=bds)
    print(res.x)

    print('plotresult')
    m,c = np.split(res.x,2)
    for i,ls in enumerate(L):
        tl = m*np.log(ls/1000.) + c
        print('model at {} : {}'.format(ls,tl))
        Xp = np.arange(X[i][:,0].min(),X[i][:,0].max()).reshape([-1,1])
        for j,noise in enumerate([-2,-3,-4,-5,-6,-7,-8]):
            Yp = flsmodel(np.hstack([Xp,noise*np.ones_like(Xp)]),tl)
            A[i].plot(Xp,Yp,'k--')
        A[i].set_xscale('log')
        A[i].set_yscale('log')
        F[i].savefig('combine_{}.png'.format(ls))
    return res.x

M = optfull()
