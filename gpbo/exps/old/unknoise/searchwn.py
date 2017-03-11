#optimize on draw from matern kernel

import gpbo
import gpbo.core.objectives as objectives
gpbo.core.debugoutput=True
import scipy as sp
import copy
import os
import time
from matplotlib import pyplot as plt
import pandas as pd
from gpbo.core.ESutils import accumulate as accumulate

run=True
plot=True
D=2
lb = [-1., -1.]
ub = [1., 1.]
s=1e-1

n=50
if run:
    ojfw, xmin, ymin = objectives.genmat52ojf(D,lb,ub)

    with open('results/searchwnopt.txt','w') as o:
        o.write('reported truemin x {} ; y {}'.format(xmin,ymin))



    def f(x,**ev):
        e = copy.deepcopy(ev)
        e['s'] = s
        if 'cheattrue' in ev.keys():
            if ev['cheattrue']:
                e['s'] = 0

        y, c, aux = ojfw(x, **e)

        return y,1,aux


    gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truegeneratedobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'),xmin=xmin)


    C=gpbo.core.config.pesfslearns(f,D,n,s,'results','eimle.csv')
    out = gpbo.search(C)

    #=gpbo.core.config.pesfsdefault(f_inplane,D,n,s,'results','matdrawfs.csv')
    #out = gpbo.search(C)

#if plot:
    #d0 = pd.read_csv('results/matdrawbs.csv')
    #plt.plot(accumulate(d0[' c']),d0[' truey at xrecc'],'r')

    #d1 = pd.read_csv('results/matdrawfs.csv')
    #plt.plot(accumulate(d1[' c']), d1[' truey at xrecc'], 'b')
    #plt.show()