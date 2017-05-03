import gpbo
import scipy as sp


D=2
n=10
s=1e-6

def f(x, **ev):
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2
    c = 1.
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c)
    return y + n, c, dict()


C=gpbo.core.config.eihypdefault(f,D,n,s,'results','eihyp.csv')
C.stoppara = {'nmax': 60}
C.stopfn = gpbo.core.optimize.nstopfn
out = gpbo.search(C)
print out