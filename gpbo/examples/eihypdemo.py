import gpbo
import scipy as sp


D=2
n=10
s=0.

def f(x, **ev):
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2
    c = 1.
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print('f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c))
    return y + n, c, dict()

def fdf(x, **ev):
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2
    dx0 = sp.sin(x[0])
    dx1 = sp.sin(x[1])
    F = sp.array([y,dx0,dx1])
    c = 1.
    n = sp.random.normal(size=3) * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print('f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c))
    return F + n, c, dict()

C=gpbo.core.config.eihypdefault(fdf,D,n,s,'results','eihyp.csv')
C.ojfchar['batchgrad']=True
C.stoppara = {'nmax': 50}
C.stopfn = gpbo.core.optimize.nstopfn
C.aqpara['mprior']= sp.array([2.,1.,2.])

out = gpbo.search(C)
print(out)