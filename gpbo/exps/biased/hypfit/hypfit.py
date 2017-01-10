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

#gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truehypobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))

#gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'worst.png'),atxa=1.)2
D=2
N=40
s=1e-9
nopts=10


if run[0]:
    for k in xrange(nopts):
        C = gpbo.core.config.pesbsdefault(f, D, N, s, 'results', 'hypfitbs{}.csv'.format(k+99))
        C.stopfn = gpbo.core.optimize.cstopfn
        C.stoppara = {'cmax': 1000}
        C.aqpara['traincfn'] = 'llogfull'
        #C.aqpara['cmax'] = C.stoppara['cmax']
        #C.stoppara['includeaq']=True
        out = gpbo.search(C)

if run[1]:
    for k in xrange(nopts):
        C = gpbo.core.config.eimledefault(f_inplane, D, N, 1e-9, 'results', 'hypfitei{}.csv'.format(k))
        C.stopfn = gpbo.core.optimize.cstopfn
        C.stoppara = {'cmax': 1000}
        C.stoppara['includeaq'] = True
        out = gpbo.search(C)



if run[2]:
    for k in xrange(nopts):
        C = gpbo.core.config.pesfsdefault(f_inplane, D, N, 1e-9, 'results', 'hypfitfs{}.csv'.format(k))
        C.stopfn = gpbo.core.optimize.cstopfn
        C.stoppara = {'cmax': 1000}
        C.stoppara['includeaq'] = True
        out = gpbo.search(C)




if plot:
    from matplotlib import pyplot as plt
    f,a = plt.subplots(1)

    d0 = [gpbo.optimize.readoptdata('results/hypfitbs{}.csv'.format(k)) for k in xrange(nopts)]
    d1 = [gpbo.optimize.readoptdata('results/hypfitfs{}.csv'.format(k)) for k in xrange(nopts)]
    d2 = [gpbo.optimize.readoptdata('results/hypfitei{}.csv'.format(k)) for k in xrange(nopts)]

    for k in xrange(nopts):

        a.plot(d0[k]['cacc'], d0[k]['trueyatxrecc'], 'b')
        print [d1[k]['cacc'][0],d1[k]['cacc'][len(d1[k]['cacc'])-1]]
        a.plot(d1[k]['cacc'], d1[k]['trueyatxrecc'], 'r')
        a.plot(d2[k]['cacc'], d2[k]['trueyatxrecc'], 'g')
    a.set_xscale('log')
    f.savefig('plots/out0.pdf')
    f, a = plt.subplots(1)

    xaxis = sp.linspace(0,1000,100)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([d0[k]['cacc'] for k in xrange(nopts)],
                                                         [d0[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)

    a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, 'b')

    low1, med1, upp1 = gpbo.core.ESutils.quartsirregular([d1[k]['cacc'] for k in xrange(nopts)],
                                                         [d1[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low1, upp1, facecolor='lightpink', edgecolor='lightpink', alpha=0.5)
    a.plot(xaxis, med1, 'r')

    low2, med2, upp2 = gpbo.core.ESutils.quartsirregular([d2[k]['cacc'] for k in xrange(nopts)],
                                                         [d2[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)

    a.fill_between(xaxis, low2, upp2, facecolor='lightgreen', edgecolor='lightgreen', alpha=0.5)
    a.plot(xaxis, med2, 'g')
    a.set_xscale('log')
    f.savefig('plots/out1.pdf')
    plt.show()

