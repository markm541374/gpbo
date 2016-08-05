#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import ESutils
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
import GPdc
import PES

nt=20
d=1
lb = sp.array([-1.]*d)
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt,d,lb,ub,GPdc.SQUEXP,sp.array([1.5,0.15]))

G = PES.makeG(X,Y,S,D,GPdc.SQUEXP,sp.array([0.,-1.]),sp.array([1.,1.]),18)
H = sp.vstack([i.hyp for i in G.kf])
f,a = plt.subplots(1)
a.plot(H[:,0],H[:,1],'r.')

np=100
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,V] = G.infer_diag(Xp,Dp)
z=5
f,a = plt.subplots(z+2)
for e,i in enumerate(sp.random.choice(range(m.shape[0]),size=z,replace=False)):
    s = sp.sqrt(V[i,:])
    a[e].fill_between(sup,sp.array(m[i,:]-2.*s).flatten(),sp.array(m[i,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
    a[e].plot(sup,m[i,:].flatten())
    a[e].plot(sp.array(X[:,0]).flatten(),Y,'g.')
    
[mp,Vp] = G.infer_diag_post(Xp,Dp)
s = sp.sqrt(Vp[0,:])
a[z].fill_between(sup,sp.array(mp[0,:]-2.*s).flatten(),sp.array(mp[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[z].plot(sup,mp[0,:].flatten())
a[z].plot(sp.array(X[:,0]).flatten(),Y,'g.')

Z=PES.drawmins(G,200,sp.array([-1.]),sp.array([1.]),SUPPORT=600,SLICELCB_PARA=1.)

a[z+1].hist(Z,bins=30)
plt.show()