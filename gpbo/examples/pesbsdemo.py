#optimize on draw from matern kernel

import gpbo
import gpbo.core.objectives as objectives
import scipy as sp
import copy
import os
import time
from matplotlib import pyplot as plt
import re


D=2
lb = [-1., -1.]
ub = [1., 1.]
s=1e-6


def f(x, **ev):
    #c = 1 - 0.5* ev['xa']
    c=45*sp.exp(-10.*ev['xa'])
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2.
    b = 0.1*ev['xa'] ** 2
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



C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','pesbsdemo.csv')
C.stoppara = {'tmax': 60 * 60*3}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
out = gpbo.search(C)


