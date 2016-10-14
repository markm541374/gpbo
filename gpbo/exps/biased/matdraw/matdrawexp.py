#optimize on draw from matern kernel

import gpbo
import gpbo.core.objectives as objectives
gpbo.core.debugoutput=True
import scipy as sp
import copy
import os
import time
from matplotlib import pyplot as plt
import re
import pandas as pd


run=False
plot=True
D=2
lb = [-1., -1.]
ub = [1., 1.]
s=1e-9

nopts=18
if run:
    for k in xrange(nopts):
        ojfw, xmin, ymin = objectives.genbiasedmat52ojf(D,lb,ub,0.15)

        with open('results/matdrawopt{}.txt'.format(k),'w') as o:
            o.write('reported truemin x {} ; y {}'.format(xmin,ymin))



        def f(x,**ev):
            y,c_,aux=ojfw(x,**ev)
            c = sp.exp(-1.*ev['xa'])
            return y,c,aux

        def f_inplane(x,**ev):
            e = copy.deepcopy(ev)
            e['xa']=0
            y, c, aux = f(x, **e)
            return y,c,aux


        gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truegeneratedobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'),xmin=xmin)

        gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'worst.png'),xmin=xmin,atxa=1.)

        C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','matdrawbs{}.csv'.format(k))
        C.stopfn = gpbo.core.optimize.cstopfn
        C.stoppara = {'cmax': 30}
        out = gpbo.search(C)


        C=gpbo.core.config.pesfsdefault(f_inplane,D,30,s,'results','matdrawfs{}.csv'.format(k))
        out = gpbo.search(C)


if plot:
    import seaborn as sns
    f,a = plt.subplots(1)
    y=sp.empty(nopts)
    r = re.compile('y (-?\d.\d+)')

    d0 = [gpbo.optimize.readoptdata('results/matdrawbs{}.csv'.format(k)) for k in xrange(nopts)]
    d1 = [gpbo.optimize.readoptdata('results/matdrawfs{}.csv'.format(k)) for k in xrange(nopts)]

    for k in xrange(nopts):
        txt = open('results/matdrawopt{}.txt'.format(k)).read()
        y[k] = float(r.findall(txt)[0])

        d0[k]['trueyatxrecc'] -=y[k]
        d1[k]['trueyatxrecc'] -=y[k]

        a.semilogy(d0[k]['cacc'], d0[k]['trueyatxrecc'], 'b')
        a.semilogy(d1[k]['cacc'], d1[k]['trueyatxrecc'], 'r')

    f, a = plt.subplots(1)

    xaxis = sp.linspace(0,30,100)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([d0[k]['cacc'] for k in xrange(nopts)],
                                                         [d0[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)

    a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, 'b')

    low1, med1, upp1 = gpbo.core.ESutils.quartsirregular([d1[k]['cacc'] for k in xrange(nopts)],
                                                         [d1[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low1, upp1, facecolor='lightpink', edgecolor='lightpink', alpha=0.5)
    a.plot(xaxis, med1, 'r')

    a.set_yscale('log')
    a.set_xlabel('accumulated cost')
    a.set_ylabel('regret')
    f.savefig('plots/out0.pdf')
    f,a = plt.subplots(1)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([d0[k]['cacc'] for k in xrange(nopts)],
                                                         [d0[k]['c'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, 'b')

    low1, med1, upp1 = gpbo.core.ESutils.quartsirregular([d1[k]['cacc'] for k in xrange(nopts)],
                                                         [d1[k]['c'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low1, upp1, facecolor='lightpink', edgecolor='lightpink', alpha=0.5)
    a.plot(xaxis, med1, 'r')
    a.set_xlabel('accumulated cost')
    a.set_ylabel('per-step cost')
    f.savefig('plots/out1.pdf')
    plt.show()

