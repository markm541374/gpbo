# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import ESutils
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt

#uniform test
class bob:
    def __init__(self,D):
        self.D=D
        return


X = ESutils.draw_support_inplane(bob(2), sp.array([-2,-1]),sp.array([0,3]),500,ESutils.SUPPORT_UNIFORM,1,0.5)
for i in xrange(X.shape[0]):
    plt.plot(X[i,0],X[i,1],'r.')
plt.axis([-5,5,-5,5])


#2d gp test

import GPdc

nt=34
X = ESutils.draw_support(bob(2), sp.array([-1.,-1.]),sp.array([1.,1.]),nt,ESutils.SUPPORT_UNIFORM)
D = [[sp.NaN]]*(nt)
hyp = sp.array([1.5,0.15,0.15])
kf = GPdc.gen_sqexp_k_d(hyp)
Kxx = GPdc.buildKsym_d(kf, X, D)
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T

lb = sp.array([-2.,-2.,-2.])
ub = sp.array([2.,2.,2.])
#MLEH =  GPdc.searchMLEhyp(X,Y,S,D,lb,ub,GPdc.SQUEXP,mx=10000)
G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, 2, sp.array([1.5, 0.15, 0.15])))
#np=180
#sup = sp.linspace(-1,1,np)
#Dp = [[sp.NaN]]*np
#Xp = sp.vstack([sp.array([i]) for i in sup])
#[m,v] = G.infer_diag(Xp,Dp)
#a0.plot(sup,m)
#sq = sp.sqrt(v)
#a0.fill_between(sup, sp.array(m-2.*sq).flatten(), sp.array(m+2.*sq).flatten(), facecolor='lightblue',edgecolor='lightblue')
plt.figure()

Z = ESutils.draw_support_inplane(G, sp.array([-1.,-1.]),sp.array([1.,1.]),400,ESutils.SUPPORT_LAPAPR,1,0.5)
R = ESutils.draw_min(G,Z,100)
for i in xrange(Z.shape[0]):
    plt.plot(Z[i,0],Z[i,1],'r.')
for i in xrange(R.shape[0]):
    plt.plot(R[i,0],R[i,1],'gx')
    
Z = ESutils.draw_support_inplane(G, sp.array([-1.,-1.]),sp.array([1.,1.]),400,ESutils.SUPPORT_LAPAPR,0,0.5)
R = ESutils.draw_min(G,Z,100)
print Z.shape
print R.shape
for i in xrange(Z.shape[0]):
    plt.plot(Z[i,0],Z[i,1],'r.')
for i in xrange(R.shape[0]):
    plt.plot(R[i,0],R[i,1],'gx')
plt.axis([-1,1,-1,1])
ng = 30
A = sp.empty([ng,ng])
for i in xrange(ng):
    for j in xrange(ng):
        A[i,j] = G.infer_LCB(sp.array([2*j/float(ng)-1.,-2*i/float(ng)+1.]),[[sp.NaN]],1.)[0,0]
plt.figure()
plt.imshow(A)


plt.show()


