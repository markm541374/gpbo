#!/usr/bin/env python2
#encoding: UTF-8

# Looking at how the min in a plane moves as I move the plane along hte 3rd axis

import scipy as sp
from scipy import linalg as spl
from matplotlib import pyplot as plt
import GPdc
import OPTutils
import ESutils
import DIRECT

#base dimension
d = 2
kindex = GPdc.MAT52
nt = 34
lb = sp.array([0.]+[-1.]*d)
ub = sp.array([5.]+[1.]*d)
Htrue = sp.array([1.4,4.]+[0.25]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt,d+1,lb,ub,kindex,Htrue, s=1e-8)
G = GPdc.GPcore(X,Y,S,D,GPdc.kernel(kindex,d+1,Htrue))

def ojfaugnn(x):
    return G.infer_m(x,[[sp.NaN]])[0,0]

def opt_ip(s):
    def dwrap(x,y):
        X = sp.hstack([[s],x])
        return (ojfaugnn(X),0)
    [xm,ym,ierror] = DIRECT.solve(dwrap,lb[1:],ub[1:], user_data=[], algmethod=1, maxf=12000, logfilename='/dev/null')
    print "DIRECT found: " +str([xm,ym,ierror])
    return xm

mintrue = opt_ip(0.)
minaug = sp.hstack([[0.],mintrue])

nm=32
spm = sp.logspace(-8,sp.log10(ub[0]),nm)
M = sp.empty([nm,d])
R = sp.empty(nm)
for i in xrange(nm):
    M[i,:] = opt_ip(spm[i])
    R[i] = spl.norm(M[i,:]-mintrue)

f,a = plt.subplots(2)
print spm
a[0].semilogy(spm,R,'bx-')
a[1].plot(M[:,0],M[:,1],'bx-')
a[1].plot(mintrue[0],mintrue[1],'ro')




np=100
u = sp.empty([d+1,np])
sup = sp.linspace(-1,1,np)
spp = sp.linspace(lb[0],ub[0],np)
for i in xrange(np):
    x=minaug.copy()
    x[0]=spp[i]
    u[0,i]=ojfaugnn(x)
    for j in xrange(d):
        x=minaug.copy()
        x[j+1]=sup[i]
        u[j+1,i] = ojfaugnn(x)
        
        

f,a = plt.subplots(d+1)
a[0].plot(spp,u[0,:].flatten())
for i in xrange(d):
    a[i+1].plot(sup,u[i+1,:].flatten())



plt.show()