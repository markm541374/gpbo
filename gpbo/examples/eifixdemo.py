import gpbo
import scipy as sp


D=4
n=10
s=1e-6

def f(x, **ev):
    y = -sp.cos(x[0]) - sp.cos(x[1]) -sp.cos(x[2]) - sp.cos(x[3]) + 4
    c = 1.
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print('f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c))
    return y + n, c, dict()


C=gpbo.core.config.eihypdefault(f,D,n,s,'results','eifix.csv')
C.stoppara = {'nmax': 200}
C.aqpara['dpara']['maxf']=2000
C.stopfn = gpbo.core.optimize.nstopfn
out = gpbo.search(C)
print(out)