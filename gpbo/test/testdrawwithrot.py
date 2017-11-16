# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import gpbo
from gpbo.core import ESutils as ESutils
import gpbo.core
from gpbo.core import debugoutput
debugoutput['drawlap']=True
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt

#uniform test
from gpbo.core import GPdc

#2d gp test


nt=15
X = ESutils.draw_support(2, sp.array([-1.,-1.]),sp.array([1.,1.]),nt,ESutils.SUPPORT_UNIFORM)
D = [[sp.NaN]]*(nt)
hyp = sp.array([1.5,0.45,0.45])
kf = GPdc.gen_sqexp_k_d(hyp)
Kxx = GPdc.buildKsym_d(kf, X, D)
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,nt)).T+sp.matrix(sps.norm.rvs(0,1e-3,nt)).T
S = sp.matrix([1e-6]*nt).T

lb = sp.array([-2.,-2.,-2.])
ub = sp.array([2.,2.,2.])

G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.SQUEXP, 2, sp.array([1.5, 0.45, 0.45])))



t  = sp.pi/12.
RotM = sp.array([[sp.cos(t),-sp.sin(t)],[sp.sin(t),sp.cos(t)]])
#Z = ESutils.draw_support(G, sp.array([-1.,-1.]),sp.array([1.,1.]),500,ESutils.SUPPORT_LAPAPR,para=np)
#R = ESutils.draw_min(G,Z,500)
f,a = plt.subplots()
Z2 = ESutils.draw_support(G, sp.array([-1.,-1.]),sp.array([1.,1.]),500,ESutils.SUPPORT_VARREJ,para=40,weighted=2,rotation=RotM)
R2 = ESutils.draw_min(G,Z2,500)

a.plot(Z2[:,0],Z2[:,1],'b.')
a.plot(R2[:,0],R2[:,1],'r.')
f.savefig('dbout/support.png')
