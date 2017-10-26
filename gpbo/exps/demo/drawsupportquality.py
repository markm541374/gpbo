import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from gpbo.core import objectives
import gpbo
import time
import pandas as pd
ki = gpbo.core.GPdc.SQUEXP
lb = np.array([-1.,-1.,-1.])
ub = np.array([1.,1.,1.])
Dim=3
def measures(n):
    print('measures with n={}'.format(n))
    f = 'results/pesfs.csv'
    names = (open(f).readline().strip('\n')+''.join([',q{}'.format(i) for i in range(5)])).replace(' ','')
    df = pd.read_csv(f,names=names.split(','),skiprows=1,engine='c')
    X = np.vstack([df['x0'].values[:n],df['x1'].values[:n],df['x2'].values[:n]]).T
    Y = df['y'].values[:n].reshape([X.shape[0],1])
    S = df['s'].values[:n].reshape([X.shape[0],1])
    D = [[np.NaN]]*n

    G = gpbo.core.PES.makeG(X,Y,S,D,ki, np.array([1.,0.,0.,0.]),np.array([1.,1.,1.,1.]),16)




    def kl(Z,G):
        m = Z.shape[0]
        c = np.zeros(m)
        for i in range(10):
            P = G.draw_post(Z,[[np.NaN]]*m,1200)
            A = np.argmin(P,axis=1)
            for a in A:
                c[a]+=1
        p = (c+1)/(m+np.sum(c))
        q = 1./float(m)
        return np.sum(-p*np.log(q/p)),np.sum(-q*np.log(p/q)),np.sum(c==0)


    m = 1000
    Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_LAPAPROT,para=8,pad_unif=False)
    laprot = kl(Z,G)

    Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_LAPAPROT,para=8,pad_unif=False,weighted=True)
    laprotw = kl(Z,G)

    Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_SLICEEI)
    ei = kl(Z,G)

    Z = gpbo.core.ESutils.draw_support(G, lb, ub, m,gpbo.core.ESutils.SUPPORT_SLICELCB,para=2)
    lcb = kl(Z,G)

    return laprot,laprotw,ei,lcb


npts = np.linspace(10,90,91).astype(int)
nv = npts.size
meas = np.zeros([4,nv])
meas2 = np.zeros([4,nv])
unuse = np.zeros([4,nv])
for i in range(nv):
    try:
        laprot,laprotw,ei,lcb = measures(npts[i])
        print(laprot,laprotw,ei,lcb)
        meas[0,i]=laprot[0]
        meas[1,i]=ei[0]
        meas[2,i]=lcb[0]
        meas[3,i]=laprotw[0]
        meas2[0,i]=laprot[1]
        meas2[1,i]=ei[1]
        meas2[2,i]=lcb[1]
        meas2[3,i]=laprotw[1]
        unuse[0,i]=laprot[2]
        unuse[1,i]=ei[2]
        unuse[2,i]=lcb[2]
        unuse[3,i]=laprotw[2]
    except:
        pass


f,a = plt.subplots(nrows=3,ncols=1)
a[0].plot(npts,meas[0,:],'bx')
a[0].plot(npts,meas[1,:],'rx')
a[0].plot(npts,meas[2,:],'gx')
a[0].plot(npts,meas[3,:],'kx')

a[1].plot(npts,unuse[0,:],'bx')
a[1].plot(npts,unuse[1,:],'rx')
a[1].plot(npts,unuse[2,:],'gx')
a[1].plot(npts,unuse[3,:],'kx')

a[2].plot(npts,meas2[0,:],'bx')
a[2].plot(npts,meas2[1,:],'rx')
a[2].plot(npts,meas2[2,:],'gx')
a[2].plot(npts,meas2[3,:],'kx')
f.savefig('figs/quality.png')