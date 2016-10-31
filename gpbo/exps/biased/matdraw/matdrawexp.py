#optimize on draw from matern kernel

import gpbo
import gpbo.core.objectives as objectives
gpbo.core.debugoutput=False
import scipy as sp
import copy
import os
import time
from matplotlib import pyplot as plt
import re
import pandas as pd

run=[True,True,True,True]
#run=[False]*4
plot=True
D=2
lb = [-1., -1.]
ub = [1., 1.]
s=1e-9

nopts=4
if any(run):
    for k in xrange(nopts):
        ojfw, xmin, ymin = objectives.genbiasedmat52ojf(D,lb,ub,0.5)

        with open('results/matdrawopt{}.txt'.format(k),'w') as o:
            o.write('reported truemin x {} ; y {}'.format(xmin,ymin))



        def f(x,**ev):
            y,c_,aux=ojfw(x,**ev)
            c = sp.exp(-3.*ev['xa'])
            return y,c,aux

        def f_inplane(x,**ev):
            e = copy.deepcopy(ev)
            e['xa']=0
            y, c, aux = f(x, **e)
            return y,c,aux

        def f_fab(x,s):
            y,c,aux = f(x,**{'xa':s,'d':[sp.NaN],'s':0.})
            return y,c


        gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truegeneratedobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'),xmin=xmin)

        gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'worst.png'),xmin=xmin,atxa=1.)
        if run[3]:

            from gpbo.exps.thirdwrap import fab as fab
            fab.runfab(lb,ub, 1., 0., f_fab, 120, 2, ofpath='results/matdrawfab{}.csv'.format(k))

        if run[0]:
            C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','matdrawbs{}.csv'.format(k))
            C.stopfn = gpbo.core.optimize.cstopfn
            C.stoppara = {'cmax': 20}
            out = gpbo.search(C)

        if run[1]:
            C=gpbo.core.config.pesfsdefault(f_inplane,D,30,s,'results','matdrawfs{}.csv'.format(k))
            C.stopfn = gpbo.core.optimize.cstopfn
            C.stoppara = {'cmax': 20}
            out = gpbo.search(C)

        if run[2]:
            C = gpbo.core.config.eimledefault(f_inplane, D, 30, s, 'results', 'matdrawei{}.csv'.format(k))
            C.stopfn = gpbo.core.optimize.cstopfn
            C.stoppara = {'cmax': 20}
            out = gpbo.search(C)

if plot:
    f,a = plt.subplots(1)
    y=sp.empty(nopts)
    r = re.compile('y (-?\d.\d+)')

    d0 = [gpbo.optimize.readoptdata('results/matdrawbs{}.csv'.format(k)) for k in xrange(nopts)]
    d1 = [gpbo.optimize.readoptdata('results/matdrawfs{}.csv'.format(k)) for k in xrange(nopts)]
    d2 = [gpbo.optimize.readoptdata('results/matdrawei{}.csv'.format(k)) for k in xrange(nopts)]
    d3 = [gpbo.optimize.readoptdata('results/matdrawfab{}.csv'.format(k)) for k in xrange(nopts)]

    for k in xrange(nopts):
        txt = open('results/matdrawopt{}.txt'.format(k)).read()
        y[k] = float(r.findall(txt)[0])

        d0[k]['trueyatxrecc'] -=y[k]
        d1[k]['trueyatxrecc'] -=y[k]
        d2[k]['trueyatxrecc'] -= y[k]
        d3[k]['trueyatxrecc'] -= y[k]

        a.semilogy(d0[k]['cacc'], d0[k]['trueyatxrecc'], 'b')
        a.semilogy(d1[k]['cacc'], d1[k]['trueyatxrecc'], 'r')
        a.semilogy(d2[k]['cacc'], d2[k]['trueyatxrecc'], 'g')
        a.semilogy(d3[k]['cacc'], d3[k]['trueyatxrecc'], 'c')
    f.savefig('plots/out0.pdf')
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

    low2, med2, upp2 = gpbo.core.ESutils.quartsirregular([d2[k]['cacc'] for k in xrange(nopts)],
                                                         [d2[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low2, upp2, facecolor='lightgreen', edgecolor='lightgreen', alpha=0.5)
    a.plot(xaxis, med2, 'g')

    low3, med3, upp3 = gpbo.core.ESutils.quartsirregular([d3[k]['cacc'] for k in xrange(nopts)],
                                                         [d3[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low3, upp3, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med3, 'c')

    a.set_yscale('log')
    a.set_xlabel('accumulated cost')
    a.set_ylabel('regret')
    f.savefig('plots/out1.pdf')
    f,a = plt.subplots(1)
    low0, med0, upp0 = gpbo.core.ESutils.quartsirregular([d0[k]['cacc'] for k in xrange(nopts)],
                                                         [d0[k]['c'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low0, upp0, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med0, 'b')

    low1, med1, upp1 = gpbo.core.ESutils.quartsirregular([d1[k]['cacc'] for k in xrange(nopts)],
                                                         [d1[k]['c'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low1, upp1, facecolor='lightpink', edgecolor='lightpink', alpha=0.5)
    a.plot(xaxis, med1, 'r')

    low2, med2, upp2 = gpbo.core.ESutils.quartsirregular([d2[k]['cacc'] for k in xrange(nopts)],
                                                         [d2[k]['c'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low2, upp2, facecolor='lightgreen', edgecolor='lightgreen', alpha=0.5)
    a.plot(xaxis, med2, 'g')

    low3, med3, upp3 = gpbo.core.ESutils.quartsirregular([d3[k]['cacc'] for k in xrange(nopts)],
                                                         [d3[k]['c'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low3, upp3, facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
    a.plot(xaxis, med3, 'c')


    a.set_xlabel('accumulated cost')
    a.set_ylabel('per-step cost')
    f.savefig('plots/out2.pdf')
    plt.show()

