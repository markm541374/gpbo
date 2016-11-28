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
n=37

def topolar(x,o):
    z = x-o
    r = sp.log10(sp.sqrt(z[0]**2+z[1]**2)) if z[0]!=0 and z[1]!=0 else -4.
    th = sp.arctan(z[1]/z[0]) if z[0] != 0. else sp.pi/2.
    return [r,1.]

ER,M,V,Z_,Y_,Ro,Y,xmin,ymin,persist = pickle.load(open("dbout/{}.p".format(n), "rb"))


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
plt.show()
