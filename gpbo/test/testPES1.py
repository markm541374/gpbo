#!/usr/bin/env python2
#encoding: UTF-8

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import ESutils
import DIRECT
import scipy as sp
from scipy import linalg as spl
from scipy import stats as sps
from matplotlib import pyplot as plt
import GPdc
import PES

nt=12
d=1
lb = sp.array([-1.]*d)
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt, d, lb, ub, GPdc.SQUEXP, sp.array([1.5, 0.15]))

G = PES.makeG(X, Y, S, D, GPdc.SQUEXP, sp.array([0., -1.]), sp.array([1., 1.]), 12)
Z=PES.drawmins(G,8,sp.array([-1.]),sp.array([1.]),SUPPORT=400,SLICELCB_PARA=1.)

Ga = GPdc.GPcore(*PES.addmins(G, X, Y, S, D, Z[0, :]) + [G.kf])

np=100
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp = sp.vstack([sp.array([i]) for i in sup])

[m,V] = G.infer_diag_post(Xp,Dp)
[mp,Vp] = Ga.infer_diag_post(Xp,Dp)

f,a = plt.subplots(2)
s = sp.sqrt(V[0,:])
a[0].fill_between(sup,sp.array(m[0,:]-2.*s).flatten(),sp.array(m[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[0].plot(sup,m[0,:].flatten())
a[0].plot(sp.array(X[:,0]).flatten(),Y,'g.')

s = sp.sqrt(Vp[0,:])
a[1].fill_between(sup,sp.array(mp[0,:]-2.*s).flatten(),sp.array(mp[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[1].plot(sup,mp[0,:].flatten())
a[1].plot(sp.array(X[:,0]).flatten(),Y,'g.')
a[1].plot(Z[0,:].flatten(),[0],'r.')

#-------------------------------------------------------------------------
#3d
nt=50
d=3
lb = sp.array([-1.]*d)
ub = sp.array([1.]*d)
[X,Y,S,D] = ESutils.gen_dataset(nt, d, lb, ub, GPdc.SQUEXP, sp.array([1.5, 0.35, 0.25, 0.30]))

G = PES.makeG(X, Y, S, D, GPdc.SQUEXP, sp.array([0., -1., -1., -1.]), sp.array([1., 1., 1., 1.]), 6)
nz=8
Z=PES.drawmins(G,nz,sp.array([-1.]*d),sp.array([1.]*d),SUPPORT=400,SLICELCB_PARA=1.)
print Z
Ga = [GPdc.GPcore(*PES.addmins(G, X, Y, S, D, Z[i, :]) + [G.kf]) for i in xrange(nz)]



np=150
sup = sp.linspace(-1,1,np)
Dp = [[sp.NaN]]*np
Xp0 = sp.vstack([sp.array([i,Z[0,1],Z[0,2]]) for i in sup])
Xp1 = sp.vstack([sp.array([Z[0,0],i,Z[0,2]]) for i in sup])
Xp2 = sp.vstack([sp.array([Z[0,0],Z[0,1],i]) for i in sup])

[mp0,Vp0] = Ga[0].infer_diag_post(Xp0,Dp)
[mp1,Vp1] = Ga[0].infer_diag_post(Xp1,Dp)
[mp2,Vp2] = Ga[0].infer_diag_post(Xp2,Dp)

[m0,V0] = G.infer_diag_post(Xp0,Dp)
[m1,V1] = G.infer_diag_post(Xp1,Dp)
[m2,V2] = G.infer_diag_post(Xp2,Dp)

f,a = plt.subplots(3)
f2,a2 = plt.subplots(3)

s = sp.sqrt(Vp0[0,:])
a[0].fill_between(sup,sp.array(mp0[0,:]-2.*s).flatten(),sp.array(mp0[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[0].plot(sup,mp0[0,:].flatten())
a[0].plot(Z[0,0].flatten(),[0],'r.')
a2[0].plot(sup,Vp0[0,:].flatten(),'g')
a2[0].plot(sup,V0[0,:].flatten(),'b')

s = sp.sqrt(Vp1[0,:])
a[1].fill_between(sup,sp.array(mp1[0,:]-2.*s).flatten(),sp.array(mp1[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[1].plot(sup,mp1[0,:].flatten())
a[1].plot(Z[0,1].flatten(),[0],'r.')

a2[1].plot(sup,Vp1[0,:].flatten(),'g')
a2[1].plot(sup,V1[0,:].flatten(),'b')


s = sp.sqrt(Vp2[0,:])
a[2].fill_between(sup,sp.array(mp2[0,:]-2.*s).flatten(),sp.array(mp2[0,:]+2.*s).flatten(),facecolor='lightblue',edgecolor='lightblue')
a[2].plot(sup,mp2[0,:].flatten())
a[2].plot(Z[0,2].flatten(),[0],'r.')

a2[2].plot(sup,Vp2[0,:].flatten(),'g')
a2[2].plot(sup,V2[0,:].flatten(),'b')

#H-----------------------------------------------------------------------
H0 = PES.PESgain(G,Ga,Z,Xp0,Dp,[1e-3]*len(Dp))
a[0].twinx().plot(sup,H0.flatten(),'r')
H1 = PES.PESgain(G,Ga,Z,Xp1,Dp,[1e-3]*len(Dp))
a[1].twinx().plot(sup,H1.flatten(),'r')
H2 = PES.PESgain(G,Ga,Z,Xp2,Dp,[1e-3]*len(Dp))
a[2].twinx().plot(sup,H2.flatten(),'r')

def cost(x,s):
    return 1.

def directwrap(Q,extra):
    x = sp.array([Q[:-1]])
    s = 10**Q[-1]
    acq = PES.PESgain(G,Ga,Z,x,[[sp.NaN]],[s])
    R = -acq/cost(x,s)
    return (R,0)

[xmin, miny, ierror] = DIRECT.solve(directwrap,sp.array([-1.,-1.,-1.,-4.]),sp.array([1.,1.,1.,0.]),user_data=[], algmethod=1, maxf=2000, logfilename='/dev/null')
print [xmin, miny, ierror]
plt.show()