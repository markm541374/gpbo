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
import GPflow as gpf
import tensorflow as tf
import tqdm
import pickle

rpath = 'fixedEI/results'
fpath = 'fixedEI/figs'
fnames=[i for i in os.listdir(rpath) if i.startswith('eihyp_3_500'.format(i)) ]

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

T = getT(rpath,fnames,250)

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
    r=200.
    r2 = 40000.
    r3 = 8000000.
    c0 = 3*y0/r2 + g0/r
    d0 = 2*y0/r3 + g0/r2
    c1 = 3*y1/r2 - g1/r
    d1 = -2*y1/r3 + g1/r2
    Y = ((d0*(x-b)+c0)*(x-b)**2 + (d1*(x-a)+c1)*(x-a)**2)
    return Y

def plotmodel(ax,theta):
    x = np.arange(250)
    y = ccmodel(x,theta)

    ax.plot(x,y,'r')
   # m0 = stats.norm.ppf(0.9,loc=0,scale=theta[0])-stats.norm.ppf(0.5,loc=0,scale=theta[0])

    #m2 = norm.ppf(0.9,loc=0,scale=theta[9])-norm.ppf(0.5,loc=0,scale=theta[9])
    #ax.plot(x[:10],y[:10]+m2,'r--')
    #ax.plot(x[:10],y[:10]-m2,'r--')
    ax.plot(x,y+y*theta[0],'r--')
    ax.plot(x,y-y*theta[0],'r--')
    return

def fitT(Y):
    theta0 = [0.1,Y[10:,0].min(),0.,Y[-1,0],0.5]
    #theta0 = [0.1,20.5,1.]#,1.]
    X = np.arange(Y.shape[0])
    def llk(theta):
        pa =  gammalogpdf(theta[0],1.,0.1)
        pa += gammalogpdf(theta[1],4.,1.)
        pa += gammalogpdf(theta[2],1.,1.)
        pa += gammalogpdf(theta[3],10.,6.)
        pa += gammalogpdf(theta[4],1.,6.)
        Yh = ccmodel(X,theta)
        #pa = gammalogpdf(theta[0],3,0.25)
        #pa += gammalogpdf(theta[1],3,1.)
        #pa += gammalogpdf(theta[2],1,1.)
        L=-pa
        #L=0

        #Yh = theta[0]*X+theta[1]
        for i in range(Y.shape[1]):
            #E = (Y[:,i]-Yh)/Yh-1

            #L-=np.sum(-0.5*((E[10:])/theta[0])**2-0.5*np.log(2*np.pi*theta[0]**2))
            L-= np.sum(normlogpdf(Y[:,i],Yh,Yh*theta[0]))
        return L
    bds = ((1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),
              (1e-9,np.Inf))#,(1e-9,np.Inf),(1e-9,np.Inf),(1e-9,np.Inf),
            # (1e-9,np.Inf),(1e-9,np.Inf))
    res = minimize(llk,theta0,method='L-BFGS-B',bounds=bds)
    return res.x

M = fitT(T)

print(M)


f,a = plt.subplots(nrows=2,ncols=1,figsize=[5,8])
for i in np.random.randint(0,T.shape[1],size=18):
    a[0].plot(T[:,i],'g-',linewidth=0.5)
    a[1].plot(np.cumsum(T[:,i]),'g-',linewidth=0.5)

plotmodel(a[0],M)
#tmean = np.mean(T,axis=1)
#tupper = tmean+2*np.sqrt(np.var(T,axis=1))
#tlower = tmean-(tupper-tmean)
#tmed = np.percentile(T,50, axis=1)
#tvar  = np.var(T,axis=1)
#a[0].plot(tmean,'b')
#a[0].plot(tupper,'b--')
#a[0].plot(tlower,'b--')
a[0].set_xlabel('Steps')
a[0].set_ylabel('overhead time (s)')

Xp = np.arange(0,250)
Yp = ccmodel(Xp,M)
a[1].plot(Xp,np.cumsum(Yp),'r')
S = np.sqrt(np.cumsum((M[0]*Yp)**2))
a[1].plot(Xp,np.cumsum(Yp)+S,'r--')
a[1].plot(Xp,np.cumsum(Yp)-S,'r--')


f.savefig(os.path.join(fpath,'overheads.png'))