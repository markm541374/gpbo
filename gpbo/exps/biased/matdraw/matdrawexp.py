#optimize on draw from matern kernel

import gpbo
import gpbo.core.objectives as objectives
gpbo.core.debugoutput=False
import scipy as sp
import copy
import os
import time
from matplotlib import pyplot as plt
import pandas as pd
from gpbo.core.ESutils import accumulate as accumulate

run=False
plot=True
D=2
lb = [-1., -1.]
ub = [1., 1.]
s=1e-9

n=40
if run:
    ojfw, xmin, ymin = objectives.genbiasedmat52ojf(D,lb,ub,0.5)

    with open('results/matdrawopt.txt','w') as o:
        o.write('reported truemin x {} ; y {}'.format(xmin,ymin))



    def f(x,**ev):
        y,c_,aux=ojfw(x,**ev)
        c = sp.exp(-2.*ev['xa'])
        return y,c,aux

    def f_inplane(x,**ev):
        e = copy.deepcopy(ev)
        e['xa']=0
        y, c, aux = f(x, **e)
        return y,c,aux


    gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'truegeneratedobjective' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'),xmin=xmin)

    gpbo.core.ESutils.plot2dFtofile(f,os.path.join('dbout', 'worst.png'),xmin=xmin,atxa=1.)

    C=gpbo.core.config.pesbsdefault(f,D,n,s,'results','matdrawbs.csv')
    out = gpbo.search(C)

    C=gpbo.core.config.pesfsdefault(f_inplane,D,n,s,'results','matdrawfs.csv')
    out = gpbo.search(C)

if plot:
    d0 = pd.read_csv('results/matdrawbs.csv')
    plt.plot(accumulate(d0[' c']),d0[' truey at xrecc'],'r')

    d1 = pd.read_csv('results/matdrawfs.csv')
    plt.plot(accumulate(d1[' c']), d1[' truey at xrecc'], 'b')
    plt.show()