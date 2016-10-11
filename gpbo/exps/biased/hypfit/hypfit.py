import scipy as sp
import numpy.random as npr
import pandas as pd
from matplotlib import pyplot as plt
import time
import gpbo
gpbo.core.debugoutput=False
import gpbo.core.GPdc as GPdc
from gpbo.core.ESutils import accumulate as accumulate
import os
import copy

run=False
plot = True

days = 148
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
f, a = plt.subplots(3)

f2, a2 = plt.subplots(1)
a = sp.hstack([a, a2])

a[0].plot(X, Y, 'g.')
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

    print "subsampling {} of {} at x={} aug={}".format(npts, n, x,ev['xa'])
    ps = npr.choice(range(n), size=npts, replace=False)
    Xd = sp.vstack([X[i] for i in ps])
    Yd = sp.vstack([Y[i] for i in ps])
    Sd = sp.vstack([S[i] for i in ps])
    Dd = [[sp.NaN]] * npts

    llk = GPdc.GP_LKonly(Xd, Yd, Sd, Dd, GPdc.kernel(GPdc.MAT52, 1, sp.array(hyp))).plk(pm, ps)
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

#gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truehypobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))

#gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'worst.png'),atxa=1.)2
D=2
N=40
s=1e-6

if run:
    C = gpbo.core.config.eimledefault(f, D, N, s, 'results', 'hypfitbs.csv')
    out = gpbo.search(C)

    C=gpbo.core.config.pesbsdefault(f,D,60,s,'results','hypfitbs.csv')
    out = gpbo.search(C)

    s=1e-6
    C=gpbo.core.config.pesfsdefault(f_inplane,D,N,s,'results','hypfitfs.csv')
    out = gpbo.search(C)

if plot:
    d0 = pd.read_csv('results/hypfitbs.csv')
    plt.plot(accumulate(d0[' c']),d0[' truey at xrecc'],'r')

    d1 = pd.read_csv('results/hypfitfs.csv')
    plt.plot(accumulate(d1[' c']), d1[' truey at xrecc'], 'b')
    plt.show()