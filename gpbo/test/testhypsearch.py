# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

from scipy import stats as sps
from scipy import linalg as spl
import scipy as sp
from matplotlib import pyplot as plt

import GPdc

ni = 100
kf = GPdc.kernel(GPdc.SQUEXP, 2, sp.array([1.3, 0.3, 0.2]))
X = sp.random.uniform(-1,1,size=[ni,2])
D = [[sp.NaN]]*ni
Kxx = GPdc.buildKsym_d(kf, X, D)
s = 1e-2
Y = spl.cholesky(Kxx,lower=True)*sp.matrix(sps.norm.rvs(0,1.,ni)).T+sp.matrix(sps.norm.rvs(0,s,ni)).T
S = sp.ones(ni)*s
print Y
MLEHYP = GPdc.searchMLEhyp(X, Y, S, D, sp.array([2., 2., 2.]), sp.array([-2., -2., -2.]), GPdc.SQUEXP)
print MLEHYP

MAPHYP = GPdc.searchMAPhyp(X, Y, S, D, sp.array([0., 0., 0.]), sp.array([1., 1., 1.]), GPdc.SQUEXP)
print MAPHYP