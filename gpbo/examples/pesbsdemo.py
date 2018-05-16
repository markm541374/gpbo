#example code to optimize using PES adapted for biased  variable fidelity evaluations

import gpbo
import scipy as sp

D=2
s=1e-6
n=100

#define a simple 2d objective in x which also varies with respect to the environmental variable
def f(x, **ev):
    #ev['xa'] is the environmental variable which varies from 0 (true objective) to 1. cost has an exponential decay away from xa=0
    c=45*sp.exp(-10.*ev['xa'])
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2.
    b = ev['xa'] ** 2
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            b = 0
            n = 0
    print('f inputs x:{} ev:{} outputs y:{} (b:{} n:{}) c:{}'.format(x, ev, y + n + b, b, n, c))
    return y + b + n, c, dict()

#arguments to generate default config are objective function, dimensionality, number of steps, noise variance, result directory and result filename
C=gpbo.core.config.pesbsdefault(f,D,n,s,'results','pesfs.csv')
#replace the default behaviour of stopping after n steps with stopping after a total evaluation and acquisition budget is exceeded
C.stoppara = {'tmax': 60 * 60 * 24}
C.stopfn = gpbo.core.optimize.totaltstopfn
#use the prediced average overhead over hte remaining optimization steps in the acquisition function
C.aqpara['overhead']='predict'
out = gpbo.search(C)
print(out)
