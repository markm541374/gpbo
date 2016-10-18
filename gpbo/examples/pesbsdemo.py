#optimize on draw from matern kernel

import gpbo
import gpbo.core.objectives as objectives
gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':True}
import scipy as sp
import copy
import os
import time
from matplotlib import pyplot as plt
import re
import pandas as pd


run=True
plot=True
D=2
lb = [-1., -1.]
ub = [1., 1.]
s=1e-9


nopts=1
if run:
    for k in xrange(nopts):
        with open('results/pesbsdemo{}.txt'.format(k),'w') as o:
            xmin=[0.,0.]
            ymin=0.
            o.write('reported truemin x {} ; y {}'.format(xmin,ymin))


        def f(x, **ev):
            #c = 1 - 0.5* ev['xa']
            c=sp.exp(-2.*ev['xa'])
            y = -sp.cos(x[0]) - sp.cos(x[1]) + 2. +0.1*s**2.
            b = ev['xa'] ** 2
            n = sp.random.normal() * sp.sqrt(s)
            if 'cheattrue' in ev.keys():
                if ev['cheattrue']:
                    b = 0
                    n = 0
            print 'f inputs x:{} ev:{} outputs y:{} (b:{} n:{}) c:{}'.format(x, ev, y + n + b, b, n, c)
            return y + b + n, c, dict()


        def f_inplane(x,**ev):
            e = copy.deepcopy(ev)
            e['xa']=0
            y, c, aux = f(x, **e)
            return y,c,aux



        C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','pesbsdemo{}.csv'.format(k))
        C.stopfn = gpbo.core.optimize.cstopfn
        C.stoppara = {'cmax': 35}
        C.aqpara['traincfn']='llog1d'
        out = gpbo.search(C)



if plot:
    import seaborn as sns
    f,a = plt.subplots(1)
    y=sp.empty(nopts)
    r = re.compile('y (-?\d.\d+)')

    d0 = [gpbo.optimize.readoptdata('results/pesbsdemo{}.csv'.format(k)) for k in xrange(nopts)]

    for k in xrange(nopts):
        txt = open('results/pesbsdemo{}.txt'.format(k)).read()
        y[k] = float(r.findall(txt)[0])

        d0[k]['trueyatxrecc'] -=y[k]

        a.semilogy(d0[k]['cacc'], d0[k]['trueyatxrecc'], 'b')

    f, a = plt.subplots(1)

    xaxis = sp.linspace(0,30,100)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([d0[k]['cacc'] for k in xrange(nopts)],
                                                         [d0[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)

    a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, 'b')


    a.set_yscale('log')
    a.set_xlabel('accumulated cost')
    a.set_ylabel('regret')
    f.savefig('plots/out0.pdf')
    f,a = plt.subplots(1)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([d0[k]['cacc'] for k in xrange(nopts)],
                                                         [d0[k]['c'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, 'b')

    a.set_xlabel('accumulated cost')
    a.set_ylabel('per-step cost')
    f.savefig('plots/out1.pdf')
    plt.show()

