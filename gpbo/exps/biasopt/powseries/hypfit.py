import scipy as sp
import numpy.random as npr
import pandas as pd

import time
import gpbo
gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':True,'taq':False}

import gpbo.core.GPdc as GPdc

import os
import copy

#run=[True,True,True]
run=[False,False,False]
plot = True


days = 120
t0 = time.clock()
df = pd.read_csv('data/DemandData_Historic-2015.csv')
t1 = time.clock()
print 'read time {0:e}'.format(t1 - t0)
N = df.shape[0]
n = min(N, days * 48)

dlb = 0.
dub = float(days)
print '{0:d} datapoints'.format(n)
X = sp.array([df.index.values[:n]]).T / 48.
Y = sp.array([df.indo.values[:n]]).T / 1000.
offs = sp.mean(Y)
Y -= offs
submode='det'
#plt.show()
S = sp.ones([n, 1]) * 1e-6
D = [[sp.NaN]] * n


pm = [1.,0.,0.]
ps = [1.,1.,1.]

def f(x, **ev):
    # print "called ojf: "+str(x)

    hyp = [10**(i*3) for i in x]

    t0 = time.clock()

    npts = max(1,int( (1-0.98*ev['xa'])*n ))
    if submode=='rand':
        print "subsampling {} of {} at x={} aug={}".format(npts, n, x,ev['xa'])
        pts = npr.choice(range(n), size=npts, replace=False)
        Xd = sp.vstack([X[i] for i in pts])
        Yd = sp.vstack([Y[i] for i in pts])
        Sd = sp.vstack([S[i] for i in pts])
        Dd = [[sp.NaN]] * npts
    elif submode=='det':
        print "first {} of {} at x={} aug={}".format(npts, n, x, ev['xa'])
        #pts = npr.choice(range(n), size=npts, replace=False)
        #print pts
        pts = map(int,sp.linspace(0,n-1,npts))
        Xd = sp.vstack([X[i] for i in pts])
        Yd = sp.vstack([Y[i] for i in pts])
        Sd = sp.vstack([S[i] for i in pts])
        Dd = [[sp.NaN]] * npts

    llk = GPdc.GP_LKonly(Xd, Yd, Sd, Dd, GPdc.kernel(GPdc.MAT52, 2, sp.array(hyp))).plk(pm, ps)
    t1 = time.clock()
    if llk < -1.:
        out = sp.log(-llk) + 1.
    else:
        out = -llk

    print "--->llk: {0} {1}    t: {2}".format(llk, out, t1 - t0)

    return out, t1-t0,dict()

def f_inplane(x,**ev):
    e = copy.deepcopy(ev)
    e['xa']=0
    y, c, aux = f(x, **e)
    return y,c,aux

