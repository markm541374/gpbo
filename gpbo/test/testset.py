# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from scipy import stats as sps
from scipy import linalg as spl
import scipy as sp
from matplotlib import pyplot as plt
import ESutils
import GPdc

nt=82
X = sp.matrix(sp.linspace(-1,1,nt)).T
D = [[sp.NaN]]*(nt)

hyp0 = sp.array([1.5,0.15])
hyp1 = sp.array([1.5,0.05])
hyp2 = sp.array([1.5,0.20])
kf = GPdc.kernel(GPdc.SQUEXP, 1, hyp0)

Kxx = GPdc.buildKsym_d(kf, X, D)

Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T

G = GPdc.GPcore(X, Y, S, D, [GPdc.kernel(GPdc.SQUEXP, 1, hyp0), GPdc.kernel(GPdc.SQUEXP, 1, hyp1), GPdc.kernel(GPdc.SQUEXP, 1, hyp2), GPdc.kernel(GPdc.SQUEXP, 1, hyp0)])
#G.printc()

np=100
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

m = G.infer_m(Xp,Dp)
f,a = plt.subplots(m.shape[0])
for i in xrange(m.shape[0]):
    a[i].plot(sup,m[i,:].flatten())
    a[i].plot(sp.array(X[:,0]).flatten(),Y,'g.')


#----------------------------------------------------
np=6
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,V] = G.infer_full(Xp,Dp)
m0 = G.infer_m(Xp,Dp)
print "xxx"
print m0
print m
print V
print "yyy"
[m,V] = G.infer_diag_post(Xp,Dp)
print m
print V
[m,V] = G.infer_full_post(Xp,Dp)
print m
print V
#--------------------------------------------
np=480
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,V] = G.infer_diag(Xp,Dp)
z=2
R = G.draw(Xp,Dp,z)
L = G.infer_LCB(Xp,Dp,2.)
E = G.infer_EI(Xp,Dp)
f,a = plt.subplots(m.shape[0])
for i in xrange(m.shape[0]):
    s = sp.sqrt(V[i,:])
    a[i].fill_between(sup,sp.array(m[i,:]-2.*s).flatten(),sp.array(m[i,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
    a[i].plot(sup,m[i,:],'b')
    a[i].plot(sp.array(X[:,0]).flatten(),Y,'g.')
    a[i].plot(sup,L[i,:].flatten(),'g')
    a[i].twinx().plot(sup,E[i,:].flatten(),'r')
    for j in xrange(z):
        a[i].plot(sup,R[z*i+j,:].flatten(),'c')
#print V
print G.llk()

nt=6
X = sp.matrix(sp.linspace(-1,1,nt)).T
D = [[sp.NaN]]*(nt)

hyp0 = sp.array([1.5,0.4])
hyp1 = sp.array([0.7,0.025])
hyp2 = sp.array([2.,0.075])
kf = GPdc.kernel(GPdc.SQUEXP, 1, sp.array([1.5, 0.15]))
Kxx = GPdc.buildKsym_d(kf, X, D)
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T
G = GPdc.GPcore(X, Y, S, D, [GPdc.kernel(GPdc.SQUEXP, 1, hyp0), GPdc.kernel(GPdc.SQUEXP, 1, hyp1), GPdc.kernel(GPdc.SQUEXP, 1, hyp2)])
m = G.infer_m(Xp,Dp)
mp = G.infer_m_post(Xp,Dp)
f,a = plt.subplots(m.shape[0]+1)
for i in xrange(m.shape[0]):
    a[i+1].plot(sup,m[i,:].flatten())
    a[i+1].plot(sp.array(X[:,0]).flatten(),Y,'g.')
a[0].plot(sp.array(X[:,0]).flatten(),Y,'g.')
a[0].plot(sup,mp[0,:].flatten())

[m,V] = G.infer_diag(Xp,Dp)
[mp,Vp] = G.infer_diag_post(Xp,Dp)
L = G.infer_LCB_post(Xp,Dp,2.)
E = G.infer_EI(Xp,Dp)
Ep = G.infer_EI_post(Xp,Dp)
z=2
W = G.draw(Xp,Dp,z)
U = G.draw_post(Xp,Dp,z)
f,a = plt.subplots(m.shape[0]+1)
for i in xrange(m.shape[0]):
    s = sp.sqrt(V[i,:])
    a[i+1].fill_between(sup,sp.array(m[i,:]-2.*s).flatten(),sp.array(m[i,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
    a[i+1].plot(sup,m[i,:],'b')
    a[i+1].plot(sp.array(X[:,0]).flatten(),Y,'g.')
    a[i+1].twinx().plot(sup,E[i,:],'r')
    for j in xrange(z):
        a[i+1].plot(sup,W[i*z+j,:],'r')
a[0].plot(sp.array(X[:,0]).flatten(),Y,'g.')
a[0].plot(sup,mp[0,:].flatten())
a[0].plot(sup,L.flatten())
a[0].twinx().plot(sup,Ep.flatten(),'r')
s = sp.sqrt(Vp[0,:])
a[0].fill_between(sup,sp.array(mp[0,:]-2.*s).flatten(),sp.array(mp[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
for j in xrange(z):
    a[0].plot(sup,U[j,:],'r')
plt.show()
