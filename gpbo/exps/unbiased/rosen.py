import gpbo
import scipy as sp
import os
import pandas as pd
from matplotlib import pyplot as plt
run = True
plot = True

D=2
n=40
s=1e-6

def f(x, **ev):
    u = 5. * x[ 0]
    v = 5. * x[ 1]
    a = 1.
    b = 100.
    y = 1e-3 * ((a - u) ** 2 + b * (v - u ** 2) ** 2)
    c = 1.
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c)
    return y + n, c, dict()



if not os.path.exists('rosen'):
    os.mkdir('rosen')

nruns = 2
if run:
    for i in range(nruns):
        C=gpbo.core.config.eimledefault(f,D,n,s,'rosen','eimle{}.csv'.format(i))
        out = gpbo.search(C)
    for i in range(nruns):
        C = gpbo.core.config.pesfsdefault(f, D, n, s, 'rosen', 'pesfs{}.csv'.format(i))
        out = gpbo.search(C)

if plot:
    f,a = plt.subplots(1)
    for i in range(nruns):
        data = pd.read_csv('rosen/eimle{}.csv'.format(i))
        a.plot( data[' truey at xrecc'],'b')

    for i in range(nruns):
        data = pd.read_csv('rosen/pesfs{}.csv'.format(i))
        a.plot( data[' truey at xrecc'],'r')

    a.set_yscale('log')

plt.show()