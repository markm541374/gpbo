#!/usr/bin/env python2
#encoding: UTF-8

#test the 1d periodic matern5/2 kernel
import sys
sys.path.append("./..")
import scipy as sp
from matplotlib import pyplot as plt
import seaborn as sns

import GPdc






kp = GPdc.kernel(GPdc.MAT52, 1, sp.array([0.75, 0.1]))
kq = GPdc.kernel(GPdc.MAT52PER, 1, sp.array([0.75, 0.1, 4.5]))
kr = GPdc.kernel(GPdc.MAT52PPT, 1, sp.array([0.75, 0.1, 4.5, 15.]))
ks = GPdc.kernel(GPdc.DEV, 1, sp.array([0.75, 0.1, 4.5, 0.5, 0.4, 10.]))
#support
ns = 1201
xax = sp.array([sp.linspace(-15,15,ns)]).T
d0 = [[sp.NaN]]
x0 = sp.array([[0.]])

p = sp.empty(ns)
q = sp.empty(ns)
r = sp.empty(ns)
s = sp.empty(ns)
for i in xrange(ns):
    p[i] = kp(x0,xax[i,:],d0,d0)
    q[i] = kq(x0,xax[i,:],d0,d0)+1
    r[i] = kr(x0,xax[i,:],d0,d0)+2
    s[i] = ks(x0,xax[i,:],d0,d0)+3
    
    
    
    
plt.plot(xax,p)
plt.plot(xax,q)
plt.plot(xax,r)
plt.plot(xax,s)


plt.show()
