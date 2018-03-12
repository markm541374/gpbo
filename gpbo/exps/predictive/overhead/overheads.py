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

def sumton2(a,n):
    #sum from i=0 to n (x-a)^2
    return (1/6.)*(n+1)*(6*a**2 - 6*a*n + 2*n**2 + n)
def sumton3(a,n):
    #sum from i=0 to n (x-a)^3
    return -0.25*(n+1)*(2*a-n)*(2*a**2 - 2*a*n + n**2 + n)
def sumton4(a,n):
    #sum from i=0 to n (x-a)^4
    return (1/30.)*(n + 1)*(30*a**4 - 60*(a**3)* n + 30. *(a**2)* n*(2* n + 1) - 30.* a* (n**2)* (n + 1) + n *(6 *n**3 + 9 *n**2 + n - 1))
def sumton5(a,n):
    #sum from i=0 to n (x-a)^5
    return  -(1/12.)* (n + 1) *(2* a - n) *(6* a**4 - 12.* (a**3)* n + 2* (a**2)* n* (7* n + 5) - 2 *a* (n**2)* (4* n + 5) + n *(2* n**3 + 4* n**2 + n - 1))
def sumton6(a,n):
    #sum from i=0 to n (x-a)^6
    a2 = a**2
    a3,a4 = a*a2, a2**2
    a5, a6 = a3*a2, a3**2
    n2 = n**2
    n3,n4 = n*n2, n2**2
    n5, n6 = n3*n2, n3**2
    return 1/42.*(n + 1)*(42.*a6 - 126.*a5*n + 210.*a4*n2 + 105.*a4*n - 210.*a3*n3 - 210.*a3*n2 + 126.*a2*n4 + 189.*a2*n3 + 21.*a2*n2 - 21.*a2*n - 42.*a*n5 - 84.*a*n4 - 21.*a*n3 + 21.*a*n2 + 6.*n6 + 15.*n5 + 6.*n4 - 6.*n3 - n2 + n)
def sumton3_3(a,b,n):
    #sum from i=0 to n (x-a)^3 (x-b)^3
    a2=a**2
    a3=a*a2
    b2=b**2
    b3=b*b2
    n2 = n**2
    n3,n4 = n*n2, n2**2
    n5,n6 = n*n4, n3**2
    return 1/420.*(n + 1)*(420.*a3*b3 - 630.*a3*b2*n + 420.*a3*b*n2 + 210.*a3*b*n - 105.*a3*n3 - 105.*a3*n2 - 630.*a2*b3*n + 1260.*a2*b2*n2 + 630.*a2*b2*n - 945.*a2*b*n3 - 945.*a2*b*n2 + 252.*a2*n4 + 378.*a2*n3 + 42.*a2*n2 - 42.*a2*n + 420.*a*b3*n2 + 210.*a*b3*n - 945.*a*b2*n3 - 945.*a*b2*n2 + 756.*a*b*n4 + 1134.*a*b*n3 + 126.*a*b*n2 - 126.*a*b*n - 210.*a*n5 - 420.*a*n4 - 105.*a*n3 + 105.*a*n2 - 105.*b3*n3 - 105.*(b**3)*(n**2) + 252.*(b**2)*(n**4) + 378.*(b**2)*(n**3) + 42.*(b**2)*(n**2) - 42.*(b**2)*n - 210.*b*n5 - 420.*b*n4 - 105.*b*n3 + 105.*b*n2 + 60.*n6 + 150.*n5 + 60.*n4 - 60.*n3 - 10.*n2 + 10.*n)
def sumton2_2(a,b,n):
    #sum from i=0 to n (x-a)^2 (x-b)^2
    a2 = a**2
    b2 = b**2
    n2 = n**2
    n3,n4 = n*n2, n2**2
    return 1/30.*(n + 1)*(30.*a2*b2 - 30.*a2*b*n + 10.*a2*n2 + 5.*a2*n - 30.*a*b2*n + 40.*a*b*n2 + 20.*a*b*n - 15.*a*n3 - 15.*a*n2 + 10.*b2*n2 + 5.*b2*n - 15.*b*n3 - 15.*b*n2 + 6.*n4 + 9.*n3 + n2 - n)
def sumton2_3(a,b,n):
    #sum from i=0 to n (x-a)^2 (x-b)^3
    a2 = a**2
    b2 = b**2
    b3,b4 = b*b2, b2**2
    n2 = n**2
    n3,n4 = n*n2, n2**2
    return 1/60.*(n + 1)*(-15.*a2*(4.*b3 - 6.*b2*n + 2.*b*n*(2.*n + 1) - n2*(n + 1)) - 2.*a*n*(-30.*b3 + 30.*b2*(2.*n + 1) - 45.*b*n*(n + 1) + 2.*(6.*n3 + 9.*n2 + n - 1)) + n*(-10.*b3*(2.*n + 1) + 45.*b2*n*(n + 1) - 6.*b*(6.*n3 + 9.*n2 + n - 1) + 5.*n*(2.*n3 + 4.*n2 + n - 1)))

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
support=390.
def ccmodel(x,theta):
    y0 = theta[1]
    g0 = theta[2]
    y1 = theta[3]
    g1 = theta[4]
    a=0.
    b=support
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
    b=support
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
    #v3d = np.einsum('ij,ij->i', np.dot(T2, v2), T2)

    v4 = v3+np.diag((m[0]*m3)**2)

    return m3,v4


def cvimodel(x,M):
    m,v = M
    a=0.
    b=support
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
    #T2 = np.hstack([((X+1-b)**3-0*(1-b)**3)/3., 0.25*((X+1-b)**4-0*(1-b)**4), ((X+1-a)**3 - 0*(1-a)**3)/3., 0.25*((X+1-a)**4-0*(1-a)**4)])
    T2 = np.hstack([sumton2(b,X), sumton3(b,X), sumton2(a,X), sumton3(a,X)])

    m3 = T2.dot(m2)
    #v3 = T2.dot(v2.dot(T2.T))
    v3d = np.einsum('ij,ij->i', np.dot(T2, v2), T2)
    coef = np.outer(m2,m2)
    #Tm = np.array([[coef[0,0]*sumton4(b,X), coef[0,1]*2*sumton5(b,X), coef[0,2]*2*sumton2_2(a,b,X), coef[0,3]*2*sumton2_3(b,a,X)],
    #               [0.          ,   coef[1,1]*sumton6(b,X), coef[1,2]*2*sumton2_3(a,b,X), coef[1,3]*2*sumton3_3(a,b,X)],
    #               [0.          , 0.            ,   coef[2,2]*sumton4(a,X)    , coef[2,3]*2*sumton5(a,X)    ],
    #               [0.          , 0.            , 0.                ,   coef[3,3]*sumton6(a,X)    ]])

    Tm = coef[0,0]*sumton4(b,X)+ coef[0,1]*2*sumton5(b,X)+ coef[0,2]*2*sumton2_2(a,b,X) + coef[0,3]*2*sumton2_3(b,a,X) + \
                                 coef[1,1]*sumton6(b,X)  + coef[1,2]*2*sumton2_3(a,b,X) + coef[1,3]*2*sumton3_3(a,b,X) + \
                                                           coef[2,2]*sumton4(a,X)       + coef[2,3]*2*sumton5(a,X)     + \
                                                                                          coef[3,3]*sumton6(a,X)


    v4 = v3d+Tm[:,0]*m[0]**2
    #print('XXXXXXXXXXXX')
    #print(Tm[-1],m[0]**2,v4[-1,-1],v3[-1,-1],v3[-1,-1]+Tm[-1]*m[0]**2)

    return m3,v4
def fitT(Y):
    assert(len(Y)>=200)
    theta0 = [0.1,Y[0],max(0.01,(Y[49]-Y[0])/50.),Y[int(support)-1],(Y[int(support)-1]-Y[0])/support]
    X = np.arange(Y.shape[0])
    def llk(theta):
        pa =  gammalogpdf(theta[0],1.,0.1)
        pa += gammalogpdf(theta[1],4.,Y[0]/4.)
        pa += gammalogpdf(theta[2],4.,max(0.01,(Y[49]-Y[0])/50.)/4.)
        pa += gammalogpdf(theta[3],4.,Y[int(support)-1]/4.)
        pa += gammalogpdf(theta[4],4.,(Y[int(support)-1]-Y[0])/support/4.)
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
    rpath = 'fixedPES/results'
    fpath = 'fixedPES/figs'
    fnames=[i for i in os.listdir(rpath) if i.startswith('pes_3_500'.format(i)) ][-8:]
    T = getT(rpath,fnames,250)
    M = buildmodel(T)
    s = 60
    B = 250*s
    plotfull(T,M,B,s)