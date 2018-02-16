import numpy as np
import scipy as sp
from scipy.optimize import minimize
from scipy.stats import gamma as gammadist
from scipy.stats import norm
from scipy.stats import lognorm
from scipy import stats
import pandas as pd
import os
import sys
from matplotlib import pyplot as plt
import tqdm
import pickle

rpath = 'fixedPES/results'
fpath = 'fixedPES/figs'
fnames=[i for i in os.listdir(rpath) if i.startswith('pes_3_500'.format(i)) ][-8:]

def getT(fpath,fnames,N):
    n = len(fnames)
    T = np.zeros([N-10,n])

    for j in tqdm.tqdm(range(n)):
        f = os.path.join(fpath,fnames[j])

        names = (open(f).readline().strip('\n')+''.join([',q{}'.format(i) for i in range(5)])).replace(' ','')
        df = pd.read_csv(f,names=names.split(','),skiprows=1,engine='c')
        tv = np.array(df['taq'].values)
        #T[:,j] = np.cumsum(tv)
        T[:,j] = tv[10:N]
    return T


gammalogpdf = lambda x,shape,scale: (shape-1)*np.log(x)-x/scale#-shape*np.log(scale)-np.log(gamma(shape))
normlogpdf = lambda x,loc,scale: -0.5*((x-loc)/scale)**2-0.5*np.log(2*np.pi*scale**2)
lognormlogpdf = lambda x,s,loc,scale: -((np.log(x)-loc)/scale)**2 - np.log(x)
tlogpdf = lambda x,v,scale: -0.5*(v+1)*np.log(1+((x/scale)**2)/v)

def ccmodel(x,theta):
    y0 = theta[1]
    g0 = theta[2]
    y1 = theta[3]
    g1 = theta[4]
    a=0.
    b=200.
    r=b-a
    r2 = r**2
    r3 = r**3
    c0 = 3*y0/r2 + g0/r
    d0 = 2*y0/r3 + g0/r2
    c1 = 3*y1/r2 - g1/r
    d1 = -2*y1/r3 + g1/r2
    Y = ((d0*(x-b)+c0)*(x-b)**2 + (d1*(x-a)+c1)*(x-a)**2)
    return Y

def cvmodel(x,M):
    m,v = M
    a=0.
    b=200.
    r=b-a
    r2=r**2
    r3=r**3
    T1 = np.array([[3./r2, 1./r,  0.,     0.],
                  [2./r3, 1./r2, 0.,     0.],
                  [0.,    0.,    3./r2, -1./r],
                  [0.,    0.,   -2./r3,  1./r2]])
    m2 = T1.dot(m[1:])
    v2 = T1.dot(v[1:,1:].dot(T1.T))

    X = x.reshape([-1,1])
    T2 = np.hstack([(X-b)**2, (X-b)**3, (X-a)**2, (X-a)**3])

    m3 = T2.dot(m2)
    v3 = T2.dot(v2.dot(T2.T))

    v4 = v3+np.diag((m[0]*m3)**2)

    return m3,v4


def fitT(Y):
    assert(len(Y)>=200)
    theta0 = [0.1,Y[0],max(0.01,(Y[49]-Y[0])/50.),Y[199],(Y[199]-Y[0])/200.]
    X = np.arange(Y.shape[0])
    def llk(theta):
        pa =  gammalogpdf(theta[0],1.,0.1)
        pa += gammalogpdf(theta[1],4.,Y[0]/4.)
        pa += gammalogpdf(theta[2],4.,max(0.01,(Y[49]-Y[0])/50.)/4.)
        pa += gammalogpdf(theta[3],4.,Y[199]/4.)
        pa += gammalogpdf(theta[4],4.,(Y[199]-Y[0])/200./4.)
        Yh = ccmodel(X,theta)
        L=0#-pa
        L-= np.sum(normlogpdf(Y,Yh,Yh*theta[0]))
        return L
    bds = ((1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),
              (1e-9,np.Inf))
    res = minimize(llk,theta0,method='L-BFGS-B',bounds=bds)
    return res.x

def buildmodel(data):
    M = [fitT(data[:,i]) for i in range(data.shape[1])]

    mean = np.mean(M,axis=0)
    print(mean)
    cov = np.cov(np.vstack(M).T)
    return [mean,cov]

def stepprobs(s,b,M):
    xmax = int(b/s)
    Xp = np.arange(xmax)
    pm,pv = cvmodel(Xp,M)
    std = np.sqrt(np.array([np.sum(pv[:i,:i]) for i in range(pv.shape[0])]))
    stepbudget = b-Xp*s
    pgx = sp.stats.norm.cdf(stepbudget,loc=np.cumsum(pm),scale=std)
    P = np.hstack([-np.diff(pgx),0])
    return P


def plotfull(T,M,B,s):
    f,a = plt.subplots(nrows=2,ncols=1,figsize=[5,7])
    for i in np.random.randint(0,T.shape[1],size=118):
        a[0].plot(T[:,i],'g-',linewidth=0.5)
        a[1].plot(np.cumsum(T[:,i]),'g-',linewidth=0.5)

    Xp = np.arange(T.shape[0])
    pm,pv = cvmodel(Xp,M)
    std = np.sqrt(np.array([np.sum(pv[:i,:i]) for i in range(pv.shape[0])]))

    a[0].plot(Xp,pm,'r')
    a[0].plot(Xp,pm+2*np.sqrt(np.diagonal(pv)),'r--')
    a[0].plot(Xp,pm-2*np.sqrt(np.diagonal(pv)),'r--')


    a[1].plot(B-s*np.arange(int(B/s)),'k')
    a[1].plot(Xp,np.cumsum(pm),'r')
    a[1].plot(Xp,np.cumsum(pm)+2*std,'r--')
    a[1].plot(Xp,np.cumsum(pm)-2*std,'r--')

    a[1].twinx().plot(stepprobs(s,B,M),'b')

    a[0].set_xlabel('Steps')
    a[0].set_ylabel('overhead time (s)')

    f.savefig(os.path.join(fpath,'overheads.png'))
    return
if __name__=="__main__":
    T = getT(rpath,fnames,250)
    M = buildmodel(T)
    s = 60
    B = 250*s
    plotfull(T,M,B,s)