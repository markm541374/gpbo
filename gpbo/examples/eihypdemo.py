import gpbo
import scipy as sp


D=4
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
f = gpbo.core.objectives.colville
q=f([0.1,0.1,0.1,0.1],**{})
C=gpbo.core.config.eihypdefault(f,D,n,s,'results','eihyp3.csv')
#C=gpbo.core.config.pesfsdefault(f,D,n,s,'results','pesfs.csv')
C.aqpara['nrandinit']=C.reccpara['onlyafter']=20
#C.ojfchar['batchgrad']=True
C.stoppara = {'nmax': 250}
C.stopfn = gpbo.core.optimize.nstopfn
#C.aqpara['mprior']= sp.array([2.,1.,2.])

out = gpbo.search(C)
print(out)