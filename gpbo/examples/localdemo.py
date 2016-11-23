import gpbo
import scipy as sp


D=2
n=100
s=0.

def f(x, **ev):
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2
    c = 1.
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c)
    return y + n, c, dict()

class conf():
    """
    fixed s, space is [-1,1]^D

    """
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.splocalaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'start':[0.65,-0.55]
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.localstopfn

        self.reccfn = gpbo.core.reccomenders.argminrecc
        self.reccpara = {
                'ev':self.aqpara['ev'],
                'check':True
                }
        self.ojfchar = {'dx': len(self.aqpara['lb']), 'dev': len(self.aqpara['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

C=conf(f,D,n,s,'results','local.csv')
out = gpbo.search(C)
print out