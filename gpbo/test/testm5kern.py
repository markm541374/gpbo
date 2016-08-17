#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
import sys
sys.path.append("./..")
import scipy as sp
from matplotlib import pyplot as plt
import seaborn as sns

import GPdc
D=4



ke = GPdc.kernel(GPdc.SQUEXP, D, sp.array([0.75, 0.1, 0.2, 0.25, 0.1]))
km = GPdc.kernel(GPdc.MAT52, D, sp.array([0.75, 0.1, 0.2, 0.25, 0.1]))
#support
ns = 1201
xax = sp.linspace(-1,1,ns)
Xo = sp.vstack([0.,0.,0.,0.]*ns)
X0 = sp.vstack([[i,0.,0.,0.] for i in xax])
X1 = sp.vstack([[0.,i,0.,0.] for i in xax])
X2 = sp.vstack([[0.,0.,i,0.] for i in xax])
X3 = sp.vstack([[0.,0.,0.,i] for i in xax])
C0 = sp.vstack([[-i,i,-i,i] for i in xax])
C1 = sp.vstack([[i,-i,i,-i] for i in xax])
C2 = sp.vstack([[i,-i,-i,i] for i in xax])
C3 = sp.vstack([[-i,i,i,-i] for i in xax])
Xax=[X0,X1,X2,X3]
Cax=[C0,C1,C2,C3]

f,a = plt.subplots(5,D)
#noderivative  --------------------------------------------------------------
d0 = [[sp.NaN]]
d1 = [[sp.NaN]]
a[0,0].set_ylabel('[n][n]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')

#firstderivative0  --------------------------------------------------------------
d0 = [[sp.NaN]]
d1 = [0]
a[1,0].set_ylabel('[n][0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#firstderivative1  --------------------------------------------------------------
d0 = [[sp.NaN]]
d1 = [1]
a[2,0].set_ylabel('[n][1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#firstderivative2  --------------------------------------------------------------
d0 = [[sp.NaN]]
d1 = [2]
a[3,0].set_ylabel('[n][2]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')
    
#firstderivative3  --------------------------------------------------------------
d0 = [[sp.NaN]]
d1 = [3]
a[4,0].set_ylabel('[n][3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[4,i].plot(xax,y,'b')
    a[4,i].plot(xax,ym,'r')

f,a = plt.subplots(4,D)
#different first derivatives--------------------------------------------------------------
#0,1  ----------------------------------------------------------------------------------
d0 = [0]
d1 = [1]
a[0,0].set_ylabel('[0][1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [2]
d1 = [0]
a[1,0].set_ylabel('[2][0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [1]
d1 = [2]
a[2,0].set_ylabel('[1][2]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [3]
d1 = [1]
a[3,0].set_ylabel('[3][1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')

f,a = plt.subplots(4,D)
#second derivative--------------------------------------------------------------
#0,1  ----------------------------------------------------------------------------------
d0 = [sp.NaN]
d1 = [0,0]
a[0,0].set_ylabel('[n][0,0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [1]
d1 = [1]
a[1,0].set_ylabel('[1][1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [2,2]
d1 = [sp.NaN]
a[2,0].set_ylabel('[2][2]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [3]
d1 = [3]
a[3,0].set_ylabel('[3][3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')
    
#third derivative--------------------------------------------------------------
#0,1  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [0]
d1 = [0,0]
a[0,0].set_ylabel('[0][0,0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [1,1]
d1 = [1]
a[1,0].set_ylabel('[1,1][1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [2,2,2]
d1 = [sp.NaN]
a[2,0].set_ylabel('[2,2,2][n]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [sp.NaN]
d1 = [3,3,3]
a[3,0].set_ylabel('[n][3,3,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')
    
#1/2 derivative--------------------------------------------------------------
#  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [1]
d1 = [2,2]
a[0,0].set_ylabel('[1][2,2]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [2,0]
d1 = [0]
a[1,0].set_ylabel('[2,0][0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [2,2,3]
d1 = [sp.NaN]
a[2,0].set_ylabel('[2,2,3][n]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [1]
d1 = [3,3]
a[3,0].set_ylabel('[1][3,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')
    
#1/1/1 derivative--------------------------------------------------------------
#  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [1]
d1 = [2,3]
a[0,0].set_ylabel('[1][2,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [2,0]
d1 = [1]
a[1,0].set_ylabel('[2,0][1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [0,2,3]
d1 = [sp.NaN]
a[2,0].set_ylabel('[0,2,3][n]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [1]
d1 = [3,0]
a[3,0].set_ylabel('[1][3,0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')
    
#1/1/1/1 derivative--------------------------------------------------------------
#  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [1,0]
d1 = [2,3]
a[0,0].set_ylabel('[1,0][2,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [0,1,2,3]
d1 = [sp.NaN]
a[1,0].set_ylabel('[1,2,3,4][n]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [sp.NaN]
d1 = [3,2,1,0]
a[2,0].set_ylabel('[n][3,2,1,0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    

#2/1/1 derivative--------------------------------------------------------------
#  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [1,1]
d1 = [2,3]
a[0,0].set_ylabel('[1,1][2,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [3,0]
d1 = [1,0]
a[1,0].set_ylabel('[3,0][1,0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [0,2,3]
d1 = [3]
a[2,0].set_ylabel('[0,2,3][3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [sp.NaN]
d1 = [1,1,2,3]
a[3,0].set_ylabel('[n][1,1,2,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')

#2/2 derivative--------------------------------------------------------------
#  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [1,1]
d1 = [2,2]
a[0,0].set_ylabel('[1,1][2,2]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [3,0]
d1 = [0,3]
a[1,0].set_ylabel('[3,0][0,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [0,2,2]
d1 = [0]
a[2,0].set_ylabel('[0,2,2][0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [sp.NaN]
d1 = [1,0,0,1]
a[3,0].set_ylabel('[n][1,0,0,1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')

#3/1 derivative--------------------------------------------------------------
#  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [1,1]
d1 = [1,2]
a[0,0].set_ylabel('[1,1][1,2]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r')
    
#2,0  --------------------------------------------------------------
d0 = [3,0]
d1 = [0,0]
a[1,0].set_ylabel('[3,0][0,0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r')
    
#1,2  --------------------------------------------------------------
d0 = [2,2,2]
d1 = [0]
a[2,0].set_ylabel('[2,2,2][0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r')
    
#3,1  --------------------------------------------------------------
d0 = [sp.NaN]
d1 = [1,0,1,1]
a[3,0].set_ylabel('[n][1,0,1,1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r')

#4 derivative--------------------------------------------------------------
#  ----------------------------------------------------------------------------------
f,a = plt.subplots(4,D)
d0 = [1,1]
d1 = [1,1]
a[0,0].set_ylabel('[1,1][1,1]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[0,i].plot(xax,y,'b')
    a[0,i].plot(xax,ym,'r.-')
    
#2,0  --------------------------------------------------------------
d0 = [0,0]
d1 = [0,0]
a[1,0].set_ylabel('[0,0][0,0]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[1,i].plot(xax,y,'b')
    a[1,i].plot(xax,ym,'r.-')
    
#1,2  --------------------------------------------------------------
d0 = [2,2,2]
d1 = [2]
a[2,0].set_ylabel('[2,2,2][2]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[2,i].plot(xax,y,'b')
    a[2,i].plot(xax,ym,'r.-')
    
#3,1  --------------------------------------------------------------
d0 = [sp.NaN]
d1 = [3,3,3,3]
a[3,0].set_ylabel('[n][3,3,3,3]')
for i in xrange(D):
    y = [ke(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    ym = [km(Xo[j,:],Cax[i][j,:],d0,d1) for j in xrange(ns)]
    a[3,i].plot(xax,y,'b')
    a[3,i].plot(xax,ym,'r.-')
plt.show()
