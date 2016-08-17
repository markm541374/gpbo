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
import DIRECT

#-------------------------------------------------------------------------
#2d
nt=60
d=2
lb = sp.array([-1.]*d)
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt, d, lb, ub, GPdc.SQUEXP, sp.array([1.5, 0.35, 0.30]))

G = PES.makeG(X, Y, S, D, GPdc.SQUEXP, sp.array([0., -1., -1.]), sp.array([1., 1., 1.]), 6)
nz=8
Z=PES.drawmins_inplane(G,nz,sp.array([-1.]*d),sp.array([1.]*d),axis=0,value = 0.,SUPPORT=500,SLICELCB_PARA=1.)

print Z
Ga = [GPdc.GPcore(*PES.addmins_inplane(G, X, Y, S, D, Z[i, :], axis=0, value=0., MINPOLICY=PES.NOMIN) + [G.kf]) for i in xrange(nz)]

np=220
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp0 = sp.vstack([sp.array([i,Z[0,1]]) for i in sup])
Xp1 = sp.vstack([sp.array([Z[0,0],i]) for i in sup])


[mp0,Vp0] = Ga[0].infer_diag_post(Xp0,Dp)
[mp1,Vp1] = Ga[0].infer_diag_post(Xp1,Dp)

f,a = plt.subplots(d)
s = sp.sqrt(Vp0[0,:])
a[0].fill_between(sup,sp.array(mp0[0,:]-2.*s).flatten(),sp.array(mp0[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[0].plot(sup,mp0[0,:].flatten())
a[0].plot(Z[0,0].flatten(),[0],'r.')

s = sp.sqrt(Vp1[0,:])
a[1].fill_between(sup,sp.array(mp1[0,:]-2.*s).flatten(),sp.array(mp1[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[1].plot(sup,mp1[0,:].flatten())
a[1].plot(Z[0,1].flatten(),[0],'r.')

H0 = PES.PESgain(G,Ga,Z,Xp0,Dp,[1e-3]*len(Dp))
a[0].twinx().plot(sup,H0.flatten(),'r')
H1 = PES.PESgain(G,Ga,Z,Xp1,Dp,[1e-3]*len(Dp))
a[1].twinx().plot(sup,H1.flatten(),'r')
##before
[mp0_,Vp0_] = G.infer_diag_post(Xp0,Dp)
[mp1_,Vp1_] = G.infer_diag_post(Xp1,Dp)


f,a = plt.subplots(d)
s = sp.sqrt(Vp0_[0,:])
a[0].fill_between(sup,sp.array(mp0_[0,:]-2.*s).flatten(),sp.array(mp0_[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[0].plot(sup,mp0_[0,:].flatten())
#a[0].plot(Z[0,0].flatten(),[0],'r.')

s = sp.sqrt(Vp1_[0,:])
a[1].fill_between(sup,sp.array(mp1_[0,:]-2.*s).flatten(),sp.array(mp1_[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[1].plot(sup,mp1_[0,:].flatten())
#a[1].plot(Z[0,1].flatten(),[0],'r.')



f,a = plt.subplots(2*d)
a[0].plot(sup,Vp0_.flatten(),'b')
a[1].plot(sup,Vp1_.flatten(),'b')

a[0].plot(sup,Vp0.flatten(),'g')
a[1].plot(sup,Vp1.flatten(),'g')

a[0].twinx().plot(sup,sp.sign((Vp0_-Vp0).flatten()),'r')
a[1].twinx().plot(sup,sp.sign((Vp1_-Vp1).flatten()),'r')

a[2].plot(sup,mp0_.flatten(),'b')
a[3].plot(sup,mp1_.flatten(),'b')

a[2].plot(sup,mp0.flatten(),'g')
a[3].plot(sup,mp1.flatten(),'g')

a[2].twinx().plot(sup,(mp0_-mp0).flatten(),'r')
a[3].twinx().plot(sup,(mp1_-mp1).flatten(),'r')


def cost(x,s):
    return 1.
axis=0
value=0
sf = 1e-4
def directwrap(Q,extra):
    x = sp.array([Q])
    s = sf
    acq = PES.PESgain(G,Ga,Z,x,[[sp.NaN]],[s])
    R = -acq/cost(x,s)
    return (R,0)

[xmin, miny, ierror] = DIRECT.solve(directwrap,sp.array([-1.,-1.]),sp.array([1.,1.]),user_data=[], algmethod=1, maxf=2000, logfilename='/dev/null')
print [xmin, miny, ierror]

plt.show()