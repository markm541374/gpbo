#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import scipy as sp
from matplotlib import pyplot as plt

n=12
def f(x):
    return 0.1*x**2-0.5*x+1.+0.1*(0.2*abs(x)+1)*sp.random.normal()
r=[]
y=[]
for i in xrange(n):
    r.append(sorted(sp.random.uniform(0,10,40)))
    y.append(map(f,r[-1]))
    plt.plot(r[-1],y[-1],'c')

import OPTutils

X,Y,lb,ub = OPTutils.mergelines(r,y)
    
plt.fill_between(X,lb,ub,facecolor='lightgreen')
plt.plot(X,Y)
plt.show()