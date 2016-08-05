#!/usr/bin/env python2
#encoding: UTF-8

# demo the derivatives for the GP lib

import scipy as sp
from scipy import linalg as spl
from matplotlib import pyplot as plt
import GPdc

f,a = plt.subplots(3)
ns = 200
sup = sp.linspace(-1,1,200)

#points are a high noise obs, a low noise obs, a derivative obs and a second derivative obs
X = sp.array([[-0.8],[-0.25],[0.25],[0.8]])
Y = sp.array([[0.3],[-0.2],[2.5],[50.]])
S = sp.array([[1e-1],[1e-6],[1e-6],[1e-6]])
D = [[sp.NaN],[sp.NaN],[0],[0,0]]

a[0].plot([-0.8,-0.8],[0.3-sp.sqrt(1e-1),0.3+sp.sqrt(1e-1)],'r')
a[0].plot([-0.8],[0.3],'ro')
a[0].plot([-0.25],[-0.2],'ro')
a[1].plot([0.25],[2.5],'ro')
a[2].plot([0.8],[50],'ro')

k= GPdc.kernel(GPdc.SQUEXP,1,sp.array([0.5,0.2]))
K = sp.empty([4,4])
for i in xrange(4):
    for j in xrange(i,4):
        K[i,j] =K[j,i] = k(X[i,:],X[j,:],D[i],D[j])
    K[i,i]+=1e-6

    
g = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,1,sp.array([0.5,0.2])))
#g.printc()
C= g.get_cho()
print C

print spl.cho_solve((C,True),sp.eye(4)).dot(K)
for i,d in enumerate([[sp.NaN],[0],[0,0]]):
    m,v = g.infer_diag(sup,[d]*ns)
    vl = (m-sp.sqrt(v)).flatten()
    vu = (m+sp.sqrt(v)).flatten()
    a[i].plot(sup,m.flatten())
    a[i].fill_between(sup,vl,vu,facecolor='LightBlue',edgecolor='LightBlue')
    a[i].set_ylabel(''.join(['d/dx ']*i)+'f')

plt.show()