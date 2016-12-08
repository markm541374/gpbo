import pickle
import scipy as sp
import numpy as np
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
from matplotlib import patches
import gpbo
from gpbo.core import GPdc as GP
from gpbo.core.choosers import gmmclassifier
from sklearn.neighbors import kneighbors_graph
from sklearn import cluster
from sklearn.preprocessing import StandardScaler
from gpbo.core import GPdc
from gpbo.core.GPdc import kernel
from gpbo.core import choosers
"""
para = {
    'lrf':lambda x:sp.exp(-2*x)+0.001*sp.exp(-0.01*x),
    'grf':lambda x:4*sp.exp(-0.4*x)+2.,
    'bcf':lambda x:0.*(7e-8)*(x+192.)**3,
    'evc':10.,
    'lsc':0.1,
    'lsn':20,
    'lsr':1e-7,
    'brm':700
}
print choosers.choice(para)

n=80

def topolar(x,o):
    z = x-o
    r = sp.log10(sp.sqrt(z[0]**2+z[1]**2)) if z[0]!=0 and z[1]!=0 else -4.
    th = sp.arctan(z[1]/z[0]) if z[0] != 0. else sp.pi/2.
    return [r,1.]
"""
n=57
ER,M,V,Z_,Y_,Ro,Y,xmin,ymin,persist = pickle.load(open("dbout/{}.p".format(n), "rb"))
Gpred = choosers.predictforward(persist['GBound'])
Lpred = choosers.predictforward(persist['LRegret'])
linit = Lpred.predict(57)
ginit = Gpred.predict(57)

    #check the switch to lacl asecision
chpara = {
    'lrf':lambda x:min(linit,Lpred.predict(x+57)),
    'grf':lambda x:min(ginit,Gpred.predict(x+57)),
    'bcf':lambda x:0.,
    'evc':1.,
    'lsc':0.,
    'lsn':20,
    'lsr':1e-7,
    'brm':90-57,
    'tol':1e-5
    }
print choosers.choice2(chpara)
"""
f,a=plt.subplots(2,2)
xaxis = sp.linspace(0,min(2*n,n+20),200)
dplot = [[sp.NaN]]*200
print persist['GBound']
for ax,data in zip([a[0,0],a[1,0],a[1,1]],[persist['ERegret'],persist['GBound'],persist['LRegret']]):

    n=len(data)
    Y = sp.array([sp.log10(i) for i in data]).reshape([n,1])
    X = sp.array([[float(i) for i in range(n)]]).T
    D = [[sp.NaN]]*n
    S = sp.array([[1e-2]]*n).T
    MAPHYP = GPdc.searchMAPhyp(X, Y, S, D, sp.array([1., 1.,0.]), sp.array([1., 1.,0.5]), GPdc.CPDEC1)
    k = kernel(GPdc.CPDEC1,1,MAPHYP)
    g = GPdc.GPcore(X, Y, S, D,k)
    m,v = g.infer_diag(xaxis,dplot)
    s = sp.sqrt(v)
    ax.fill_between(xaxis,(m-2.*s).flatten(),(m+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
    ax.plot(xaxis,m.flatten(),'b')
    ax.plot(X,Y,'r.')

    #ax.set_yscale('log')
ns = Ro.shape[0]
R = sp.empty(Ro.shape)
for i in xrange(ns):
    R[i,:] = topolar(Ro[i,:],xmin)
#    print Ro[i,:],R[i,:]
aug = sp.empty([R.shape[0],1])
for i in xrange(ns):
    aug[i,0] = (R[i,0]-xmin[0])**2 + (R[i,1]-xmin[1])**2


Y_,classaux = gmmclassifier(R,xmin)

F = StandardScaler().fit_transform(R)
bandwidth = cluster.estimate_bandwidth(F, quantile=0.3)
ms = cluster.MeanShift(bandwidth=bandwidth, bin_seeding=True)
ms.fit(R)
Ys_ = ms.predict(R)

f,a=plt.subplots(4)

for i in xrange(ns):
    sym = ['r.','b.', 'g.','c.','k.','bx','rx','gx','cx','kx']+['grey']*20
    col = sym[Ys_[i]]
    a[0].plot(R[i,0],R[i,1],col)
    col = sym[Ys_[i]]
    a[1].plot(Ro[i, 0], Ro[i, 1], col)

for i in xrange(ns):
    sym = ['r.','b.', 'g.','c.','k.','bx','rx','gx','cx','kx']+['grey']*20
    col = sym[Y_[i]]
    a[2].plot(R[i,0],R[i,1],col)
    col = sym[Y_[i]]
    a[3].plot(Ro[i, 0], Ro[i, 1], col)
"""
plt.show()
