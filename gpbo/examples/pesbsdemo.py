import gpbo
import scipy as sp
from gpbo.core import GPdc as GPdc
gpbo.core.debugoutput=True
gpbo.core.debugoptions={'datavis':False,'drawlap':False,'cost1d':True}
D=2
n=30
s=1e-9

def f(x,**ev):

    c=sp.exp(-2.*ev['xa'])

    y=-sp.cos(x[0])-sp.cos(x[1])+2
    b=ev['xa']**2
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            b=0
            n=0
    print 'f inputs x:{} ev:{} outputs y:{} (b:{} n:{}) c:{}'.format(x, ev, y+n+b, b, n, c)
    return y+b+n,c,dict()



C=gpbo.core.config.pesbsdefault(f,D,n,s,'results','pesbs.csv')
C.aqpara['SUPPORT_MODE']=[gpbo.core.ESutils.SUPPORT_LAPAPROT]
out = gpbo.search(C)
print out