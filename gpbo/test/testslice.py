#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from matplotlib import pyplot as plt
import slice
import scipy as sp
def llk(x):
    return sp.log( sp.exp(-(0.5*(x-1)**2)/0.5**2)+sp.exp(-(0.5*(x+2)**2)/0.25**2))

S = slice.slice_sample(llk,sp.array([0.]),5000,sp.array([0.05]))
#print S
plt.hist(sp.array(S).flatten(),bins=30)

def llk2(x):
    return sp.log( sp.exp(-0.5*((x[0]+1.5)**2+(x[1]-1.5)**2)/0.5**2)+sp.exp(-0.5*(((x[0]-1.5)**2)/0.75**2+((x[1]+1.5)**2)/0.25**2)))

S = slice.slice_sample(llk2,sp.array([0.,0.]),10000,sp.array([0.15,0.15]))
A = sp.zeros([20,20])
for i in xrange(S.shape[0]):
    a1 = int((S[i,0]+3)*20/6)
    a2 = int((S[i,1]+3)*20/6)
    try:
        A[a1,a2]+=1
    except:
        pass
#plt.figure()
#for i in xrange(A.shape[0]):
#    plt.plot(S[i,0],S[i,1],'r.')
plt.figure()
plt.imshow(A)
plt.show()