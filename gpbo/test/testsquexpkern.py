# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from scipy import stats as sps
from scipy import linalg as spl
import scipy as sp
from matplotlib import pyplot as plt

from gpbo.core import GPdc

nt=411
X = sp.matrix(sp.linspace(-1,1,nt)).T
D = [[sp.NaN]]*(nt)

hyp = sp.array([1.5,0.15])
kf = GPdc.gen_sqexp_k_d(hyp)
print "X"
#print kf(sp.array([0.1]),sp.array([0.2]),[[sp.NaN]],[[sp.NaN]])
Kxx = GPdc.buildKsym_d(kf, X, D)
print "X"
Y = spl.cholesky(Kxx+sp.eye(nt)*1e-6,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-4]*nt).T
f0,a0 = plt.subplots(1)
a0.plot(sp.array(X[:,0]).flatten(),Y,'g.')


lb = sp.array([-2.,-2.])
ub = sp.array([2.,2.])
import time
t0=time.time()
MLEH =  GPdc.searchMLEhyp(X, Y, S, D, lb, ub, GPdc.SQUEXP, mx=10000)
print "X"
print MLEH
mp = sp.array([0.,-1.])
sb = sp.array([1.,1.])
#MAPH =  GPdc.searchMAPhyp(X,Y,S,D,mp,sb,GPdc.SQUEXP,mx=10000)
G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, 1, MLEH))
print 'time {}'.format(t0-time.time())
print "X"
print Y.shape
#G.printc()
print G.llk()

np=180
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,v] = G.infer_diag(Xp,Dp)
a0.plot(sup,m.flatten())
sq = sp.sqrt(v)

a0.fill_between(sup, sp.array(m-1.*sq).flatten(), sp.array(m+1.*sq).flatten(), facecolor='lightblue',edgecolor='lightblue')
#f0.savefig('../figs/testsquexpkern.png')
plt.show()
