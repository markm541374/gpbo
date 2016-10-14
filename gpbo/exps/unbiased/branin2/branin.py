import gpbo
import scipy as sp
import os
import pandas as pd
from matplotlib import pyplot as plt
run = [True,True,True]
startfrom=0
plot = True

D=2

lsl=-7
lsu=3
nopts=1
cmax=50
def cfn(x,s):
    return s**-0.5

def f(x, **ev):
    u = x[0] * 7.5 + 2.5
    v = x[1] * 7.5 + 2.5
    y = (-1.275 * (u / sp.pi) ** 2 + 5 * u / sp.pi + v - 6) ** 2 + (10. - 5. / (4 * sp.pi)) * sp.cos(u) + 10.

    c = cfn(x, ev['s'])
    n = sp.random.normal() * sp.sqrt(ev['s'])
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
            c=42
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c)
    return y + n, c, dict()

if run[0]:
    for k in xrange(nopts):
        C=gpbo.core.config.pesvsdefault(f,cfn,D,50,lsl,lsu,'results','pesvs{}.csv'.format(k+startfrom))
        #C.stopfn = gpbo.core.optimize.cstopfn
        #C.stoppara = {'cmax': cmax}
        out = gpbo.search(C)
if run[1]:
    for k in xrange(nopts):
        C = gpbo.core.config.pesfslearns(f, D, 1, 10**lsl, 'results', 'pesfslow{}.csv'.format(k+startfrom))
        C.stopfn = gpbo.core.optimize.cstopfn
        C.stoppara = {'cmax': cmax}
        out = gpbo.search(C)
if run[2]:
    for k in xrange(nopts):
        C = gpbo.core.config.pesfslearns(f, D, 50, 10**lsu, 'results', 'pesfshigh{}.csv'.format(k+startfrom))
        #C.stopfn = gpbo.core.optimize.cstopfn
        #C.stoppara = {'cmax': cmax}
        out = gpbo.search(C)
if plot:
    f, a = plt.subplots(1)
    y = 0.39788735772973816


    d0 = [gpbo.optimize.readoptdata('results/pesvs{}.csv'.format(k)) for k in xrange(nopts)]
    d1 = [gpbo.optimize.readoptdata('results/pesfslow{}.csv'.format(k)) for k in xrange(nopts)]
    d2 = [gpbo.optimize.readoptdata('results/pesfshigh{}.csv'.format(k)) for k in xrange(nopts)]

    for k in xrange(nopts):
        d0[k]['trueyatxrecc'] -= y
        a.semilogy(d0[k]['cacc'], d0[k]['trueyatxrecc'], 'b')

        d1[k]['trueyatxrecc'] -= y
        a.semilogy(d1[k]['cacc'], d1[k]['trueyatxrecc'], 'r')

        d2[k]['trueyatxrecc'] -= y
        a.semilogy(d2[k]['cacc'], d2[k]['trueyatxrecc'], 'g')

    f, a = plt.subplots(1)
    xaxis = sp.linspace(0, cmax, 100)
    low0,med0,upp0 = gpbo.core.ESutils.quartsirregular([d0[k]['cacc'] for k in xrange(nopts)],
                                             [d0[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)

    a.fill_between(xaxis,low0,upp0,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
    a.plot(xaxis, med0, 'b')



    xaxis = sp.linspace(0, cmax, 100)
    low1,med1,upp1 = gpbo.core.ESutils.quartsirregular([d1[k]['cacc'] for k in xrange(nopts)],
                                             [d1[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low1, upp1, facecolor='lightpink', edgecolor='lightpink',alpha=0.5)
    a.plot(xaxis, med1, 'r')

    xaxis = sp.linspace(0, cmax, 100)
    low2,med2,upp2 = gpbo.core.ESutils.quartsirregular([d2[k]['cacc'] for k in xrange(nopts)],
                                             [d2[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)
    a.fill_between(xaxis, low2, upp2, facecolor='lightgreen', edgecolor='lightgreen',alpha=0.5)
    a.plot(xaxis, med2, 'g')
    a.set_yscale('log')
    plt.show()