# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import gpbo
from gpbo.core import ESutils as ESutils
import gpbo.core
gpbo.core.debugoutput=True
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt

#uniform test
from gpbo.core import GPdc

#2d gp test


nt=34
X = ESutils.draw_support(2, sp.array([-1.,-1.]),sp.array([1.,1.]),nt,ESutils.SUPPORT_UNIFORM)
D = [[sp.NaN]]*(nt)
hyp = sp.array([1.5,0.25,0.25])
kf = GPdc.gen_sqexp_k_d(hyp)
Kxx = GPdc.buildKsym_d(kf, X, D)
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T

lb = sp.array([-2.,-2.,-2.])
ub = sp.array([2.,2.,2.])

G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, 2, sp.array([1.5, 0.25, 0.25])))



np=40

Z = ESutils.draw_support(G, sp.array([-1.,-1.]),sp.array([1.,1.]),500,ESutils.SUPPORT_LAPAPR,para=np)
R = ESutils.draw_min(G,Z,500)

Z2 = ESutils.draw_support(G, sp.array([-1.,-1.]),sp.array([1.,1.]),500,ESutils.SUPPORT_LAPAPROT,para=np)
R2 = ESutils.draw_min(G,Z2,500)

plt.figure()
for i in xrange(Z.shape[0]):
    plt.plot(Z[i,0],Z[i,1],'r.')
for i in xrange(Z2.shape[0]):
    plt.plot(Z2[i,0],Z2[i,1],'b.')

plt.axis([-1,1,-1,1])
ng = 30
A = sp.empty([ng,ng])
for i in xrange(ng):
    for j in xrange(ng):
        A[i,j] = G.infer_m_post(sp.array([2*j/float(ng)-1.,-2*i/float(ng)+1.]),[[sp.NaN]])[0,0]
plt.figure()
plt.imshow(A)




plt.show()