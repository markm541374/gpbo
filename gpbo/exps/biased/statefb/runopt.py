import scipy as sp
import numpy.random as npr
import pandas as pd

import time
import gpbo
gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':False,'taq':False,'support':False}

import gpbo.core.GPdc as GPdc

import os
import copy

run=[True,True,True,True]
#run=[False,False,False,True]
plot = True

import statefb
f=statefb.f
f_inplane = statefb.f_inplane



D=2
N=40
s=1e-9
nopts=8


if run[0]:
    for k in xrange(nopts):
        C = gpbo.core.config.pesbsdefault(f, D, N, s, 'results', 'statefbbs{}.csv'.format(k))
        C.stoppara = {'tmax': 60 * 60 *2}
        C.stopfn = gpbo.core.optimize.totaltstopfn
        C.aqpara['overhead'] = 'last'
        C.aqpara['traincfn'] = 'llogfull'
        #C.aqpara['cmax'] = C.stoppara['cmax']
        out = gpbo.search(C)

if run[1]:
    for k in xrange(nopts):
        C = gpbo.core.config.eimledefault(f_inplane, D, N, 1e-9, 'results', 'statefbei{}.csv'.format(k))
        C.stoppara = {'tmax': 60 * 60*2}
        C.stopfn = gpbo.core.optimize.totaltstopfn
        C.aqpara['overhead'] = 'last'
        out = gpbo.search(C)



if run[2]:
    for k in xrange(nopts):
        C = gpbo.core.config.pesfsdefault(f_inplane, D, N, 1e-9, 'results', 'statefbfs{}.csv'.format(k))
        C.stoppara = {'tmax': 60 * 60*2}
        C.stopfn = gpbo.core.optimize.totaltstopfn
        C.aqpara['overhead'] = 'last'
        out = gpbo.search(C)
from gpbo.exps.thirdwrap.mtbowrap import optmtbo
if run[3]:
    for k in xrange(nopts):
        out = optmtbo(f,sp.array([-1.]*3),sp.array([1.]*3),0.5,50,fname='statefbmt{}.csv'.format(k),fpath='results')



if plot:
    from matplotlib import pyplot as plt
    f,a = plt.subplots(1)

    d0 = [gpbo.optimize.readoptdata('results/statefbbs{}.csv'.format(k),includetaq=True) for k in xrange(nopts)]
    d1 = [gpbo.optimize.readoptdata('results/statefbfs{}.csv'.format(k),includetaq=True) for k in xrange(nopts)]
    d2 = [gpbo.optimize.readoptdata('results/statefbei{}.csv'.format(k),includetaq=True) for k in xrange(nopts)]

    d3 = [gpbo.optimize.readoptdata('results/statefbmt{}.csv'.format(k),includetaq=True) for k in xrange(nopts)]
    for k in xrange(nopts):

        a.plot(d0[k]['cacc'], d0[k]['trueyatxrecc'], 'b')
        print [d1[k]['cacc'][0],d1[k]['cacc'][len(d1[k]['cacc'])-1]]
        a.plot(d1[k]['cacc'], d1[k]['trueyatxrecc'], 'r')
        a.plot(d2[k]['cacc'], d2[k]['trueyatxrecc'], 'g')
        a.plot(d3[k]['cacc'], d3[k]['trueyatxrecc'], 'k')
    a.set_xscale('log')
    f.savefig('plots/out0.pdf')
    f, a = plt.subplots(1)

    xaxis = sp.linspace(0,60*60*2,1000)
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

    low3, med3, upp3 = gpbo.core.ESutils.quartsirregular([d3[k]['cacc'] for k in xrange(nopts)],
                                                         [d3[k]['trueyatxrecc'] for k in xrange(nopts)], xaxis)

    a.fill_between(xaxis, low3, upp3, facecolor='lightgreen', edgecolor='lightgreen', alpha=0.5)
    a.plot(xaxis, med3, 'k')
    a.set_xscale('log')
    f.savefig('plots/out1.pdf')
    plt.show()

