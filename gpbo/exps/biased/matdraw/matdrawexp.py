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


run=True
plot=True
D=2
lb = [-1., -1.]
ub = [1., 1.]
s=1e-9

nopts=1
if run:
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


        gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truegeneratedobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'),xmin=xmin)

        gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'worst.png'),xmin=xmin,atxa=1.)

        C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','matdrawbs{}.csv'.format(k))
        C.stopfn = gpbo.core.optimize.cstopfn
        C.stoppara = {'cmax': 20}
        out = gpbo.search(C)


        C=gpbo.core.config.pesfsdefault(f_inplane,D,20,s,'results','matdrawfs{}.csv'.format(k))
        out = gpbo.search(C)


if plot:
    f,a = plt.subplots(1)
    y=sp.empty(nopts)
    r = re.compile('y (-?\d.\d+)')

    d0 = [gpbo.optimize.readoptdata('results/matdrawbs{}.csv'.format(k)) for k in xrange(nopts)]
    d1 = [gpbo.optimize.readoptdata('results/matdrawfs{}.csv'.format(k)) for k in xrange(nopts)]

    for k in xrange(nopts):
        txt = open('results/matdrawopt{}.txt'.format(k)).read()
        print txt
        y[k] = float(r.findall(txt)[0])
        d0[k]['trueyatxrecc'] -=y[k]
        d1[k]['trueyatxrecc'] -=y[k]

        a.semilogy(d1[k]['cacc'], d1[k]['trueyatxrecc'], 'b')
        a.semilogy(d0[k]['cacc'], d0[k]['trueyatxrecc'], 'r')

    xaxis = sp.linspace(0,20,100)
    med0 = gpbo.core.ESutils.medianirregular([d0[k]['cacc'] for k in xrange(nopts)],[d0[k]['trueyatxrecc'] for k in xrange(nopts)],xaxis)
    a.plot(xaxis,med0,'r.')

    med1 = gpbo.core.ESutils.medianirregular([d1[k]['cacc'] for k in xrange(nopts)], [d1[k]['trueyatxrecc'] for k in xrange(nopts)],
                           xaxis)
    a.plot(xaxis, med1, 'b.')
    plt.show()

