#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import eprop
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt

#basic ep test on three values
m = sp.array([1.,2.,0.5])
v = sp.array([[1.,0.5,0.2],[0.5,2.,0.5],[0.2,0.5,1.]])
n=len(m)


Y = sp.array([1.,0.,1.])
Z = sp.array([-1,1,-1])
F = sp.array([0.25**2,0.,0.])
mt,vt = eprop.expectation_prop(m,v,Y,Z,F,35)

mu,vu = eprop.gaussian_fusion(m,mt,v,vt)

f,ax = plt.subplots(n,sharex=True)
ns=200
sup = sp.linspace(-6,6,ns)
for i in xrange(n):
    ax[i].plot(sup,sps.norm.pdf(sup,loc=m[i],scale=sp.sqrt(v[i,i])),'b')
    ax[i].plot(sup,sps.norm.pdf(sup,loc=mu[i],scale=sp.sqrt(vu[i,i])),'r')
    ax[i].twinx().plot(sup,sps.norm.cdf(Z[i]*(sup-Y[i])/max(F[i],1e-20))*sps.norm.pdf(sup,loc=m[i],scale=sp.sqrt(v[i,i])),'g')
#est on a gp

import ESutils
import GPdc
nt=5
X = sp.matrix(sp.linspace(-1,1,nt)).T
D = [[sp.NaN]]*(nt)
hyp = sp.array([1.5,0.15])
kf = GPdc.kernel(GPdc.SQUEXP, 1, hyp)
Kxx = GPdc.buildKsym_d(kf, X, D)
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-2]*nt).T
g = GPdc.GPcore(X, Y, S, D, kf)
f,a = plt.subplots(2)

Xe = sp.array([-0.25,0.25])
De = [[sp.NaN]]*2
[m0,V0] = g.infer_full(Xe,De)
Ye = sp.array([2.,-2.])
Ze = sp.array([1.,-1.])
Fe = sp.array([(1e-6)**2,(1e-6)**2])

mt,vt = eprop.expectation_prop(m0,V0,Ye,Ze,Fe,20)
print D+De
g2 = GPdc.GPcore(sp.vstack([X, sp.array([Xe]).T]), sp.vstack([Y, sp.array([Ye]).T]), sp.vstack([S, sp.array([Fe]).T]), D + De, kf)

ESutils.plot_gp(g,a[0],sp.linspace(-1,1,100),[[sp.NaN]]*100)
ESutils.plot_gp(g2,a[1],sp.linspace(-1,1,100),[[sp.NaN]]*100)
plt.show()