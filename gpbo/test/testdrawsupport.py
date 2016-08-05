# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import ESutils
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt

#uniform test
import GPdc

X = ESutils.draw_support(2, sp.array([-2,-1]),sp.array([0,3]),500,ESutils.SUPPORT_UNIFORM)
for i in xrange(X.shape[0]):
    plt.plot(X[i,0],X[i,1],'r.')
plt.axis([-5,5,-5,5])


#2d gp test



nt=34
X = ESutils.draw_support(2, sp.array([-1.,-1.]),sp.array([1.,1.]),nt,ESutils.SUPPORT_UNIFORM)
D = [[sp.NaN]]*(nt)
hyp = sp.array([1.5,0.25,0.25])
kf = GPdc.gen_sqexp_k_d(hyp)
Kxx = GPdc.buildKsym_d(kf,X,D)
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T

lb = sp.array([-2.,-2.,-2.])
ub = sp.array([2.,2.,2.])
#MLEH =  GPdc.searchMLEhyp(X,Y,S,D,lb,ub,GPdc.SQUEXP,mx=10000)
G = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,2,sp.array([1.5,0.15,0.15])))
#np=180
#sup = sp.linspace(-1,1,np)
#Dp = [[sp.NaN]]*np
#Xp = sp.vstack([sp.array([i]) for i in sup])
#[m,v] = G.infer_diag(Xp,Dp)
#a0.plot(sup,m)
#sq = sp.sqrt(v)
#a0.fill_between(sup, sp.array(m-2.*sq).flatten(), sp.array(m+2.*sq).flatten(), facecolor='lightblue',edgecolor='lightblue')

Z = ESutils.draw_support(G, sp.array([-1.,-1.]),sp.array([1.,1.]),500,ESutils.SUPPORT_SLICEEI)
np=30
R = ESutils.draw_min(G,Z,500)
Z2 = ESutils.draw_support(G, sp.array([-1.,-1.]),sp.array([1.,1.]),500,ESutils.SUPPORT_LAPAPR,para=np)
R2 = ESutils.draw_min(G,Z2,500)

plt.figure()
for i in xrange(Z.shape[0]):
    plt.plot(Z[i,0],Z[i,1],'r.')
for i in xrange(Z2.shape[0]):
    plt.plot(Z2[i,0],Z2[i,1],'b.')
for i in xrange(R.shape[0]):
    plt.plot(R[i,0],R[i,1],'gx')
plt.axis([-1,1,-1,1])
ng = 30
A = sp.empty([ng,ng])
for i in xrange(ng):
    for j in xrange(ng):
        A[i,j] = G.infer_m_post(sp.array([2*j/float(ng)-1.,-2*i/float(ng)+1.]),[[sp.NaN]])[0,0]
plt.figure()
plt.imshow(A)


#1d gp test
nt=10
X = sp.matrix(sp.linspace(-1,1,nt)).T
D = [[sp.NaN]]*(nt)
hyp = sp.array([1.5,0.15])
kf = GPdc.gen_sqexp_k_d(hyp)
Kxx = GPdc.buildKsym_d(kf,X,D)
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T
f0 = plt.figure()
a0 = plt.subplot(311)
a1 = plt.subplot(312)
a2 = plt.subplot(313)
a0.plot(sp.array(X[:,0]).flatten(),Y,'g.')
lb = sp.array([-2.,-2.])
ub = sp.array([2.,2.])
#MLEH =  GPdc.searchMLEhyp(X,Y,S,D,lb,ub,GPdc.SQUEXP,mx=10000)
G = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,1,sp.array([1.5,0.15])))
np=180
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])
[m,v] = G.infer_diag(Xp,Dp)
a0.plot(sup,m.flatten())
sq = sp.sqrt(v)
a0.fill_between(sup, sp.array(m-2.*sq).flatten(), sp.array(m+2.*sq).flatten(), facecolor='lightblue',edgecolor='lightblue')


X = ESutils.draw_support(G, sp.array([-1]),sp.array([1.]),300,ESutils.SUPPORT_SLICEEI)
a1.hist(X,bins=80)

R = ESutils.draw_min(G,X,300)

j = a2.hist(R,bins=60)
a2.axis([-1,1,0,1.2*max(j[0])])

plt.show()