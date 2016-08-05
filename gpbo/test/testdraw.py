# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


import ctypes as ct
import scipy as sp
from scipy import stats as sps
from scipy import linalg as spl
from matplotlib import pyplot as plt
import os
libGP = ct.cdll.LoadLibrary(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../../../dist/Release/GNU-Linux/libGPshared.so'))

ctpd = ct.POINTER(ct.c_double)

import GPdc

ni = 50
kf = GPdc.kernel(GPdc.SQUEXP,1,sp.array([1.3,0.3]))
X = sp.matrix(sp.linspace(-1,1,ni)).T
D = [[sp.NaN]]*ni
Kxx = GPdc.buildKsym_d(kf,X,D)

tmp = spl.cholesky(Kxx,lower=True)
Ch = sp.vstack([tmp[i,:] for i in xrange(ni)]) #hack/force row major storage


z = 5

b = sp.empty([ni,z])


libGP.drawcov(Ch.ctypes.data_as(ctpd),ct.c_int(ni),b.ctypes.data_as(ctpd),ct.c_int(z))

#print b
for i in xrange(z):
    plt.plot(X[:,0],b[:,i])




#copy of test squexp from here
import GPdc

nt=12
X = sp.matrix(sp.linspace(-1,1,nt)).T
D = [[sp.NaN]]*(nt)

hyp = sp.array([1.5,0.15])
kf = GPdc.gen_sqexp_k_d(hyp)

Kxx = GPdc.buildKsym_d(kf,X,D)

Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T
f0 = plt.figure()
a0 = plt.subplot(111)
a0.plot(sp.array(X[:,0]).flatten(),Y,'g.')


lb = sp.array([-2.,-2.])
ub = sp.array([2.,2.])
MLEH =  GPdc.searchMLEhyp(X,Y,S,D,lb,ub,GPdc.SQUEXP,mx=10000)

print MLEH
G = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,1,MLEH))
print G.llk()

np=180
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,v] = G.infer_diag(Xp,Dp)

a0.plot(sup,m.flatten())
sq = sp.sqrt(v)

a0.fill_between(sup, sp.array(m-2.*sq).flatten(), sp.array(m+2.*sq).flatten(), facecolor='lightblue',edgecolor='lightblue')
#to here

z=5
R = G.draw(Xp,Dp,z)
print R
for i in xrange(z):
    plt.plot(sup,sp.array(R[i,:]).flatten(),'r')

plt.show()
