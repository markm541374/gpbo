#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from matplotlib import pyplot as plt
from gpbo.core.slice import slice_sample
import numpy as np
import scipy as sp
def llk(x):
    return sp.log( .5*sp.exp(-(0.5*(x-1)**2)/0.75**2)*(0.5*x+4) +sp.exp(-(0.5*(x+1)**2)/0.25**2))

S = slice_sample(llk,sp.array([0.]),10000,sp.array([0.05]))
#print S
f,a = plt.subplots()
a.hist(sp.array(S).flatten(),bins=70,normed=1)
x = np.linspace(-3,4,200)
y = np.exp(llk(x))
y/= np.sum(y)*(7/200.)
a.plot(x,y,'r')
f.savefig('slice.png')
plt.close()
