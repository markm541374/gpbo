import gpbo
import scipy as sp
from gpbo.core import GPdc

D=2
n=30
s=1e-3

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
        aq0 = gpbo.core.acquisitions.EIMAPaq
        aq0para= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D,
            'nrandinit': 10,
            'mprior': sp.array([1.]+[0.]*D),
            'sprior': sp.array([1.]*(D+1)),
            'kindex': GPdc.MAT52,
            'volper':1e-6
        }

        aq1 = gpbo.core.acquisitions.bruteaq
        aq1para = {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.] * D,
            'ub': [1.] * D,
        }

        self.aqfn = gpbo.core.acquisitions.choiceaq
        self.aqpara = {
            'aqoptions':[[aq0,aq0para],[aq1,aq1para]],
            'chooser':gpbo.core.choosers.alternate,
            'choosepara':dict(),
            'ev': aq0para['ev'],
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.nstopfn

        self.reccfn = gpbo.core.reccomenders.gpmaprecc
        self.reccpara = {
                'ev':aq0para['ev'],
                'lb':aq0para['lb'],
                'ub':aq0para['ub'],
                'mprior':aq0para['mprior'],
                'sprior':aq0para['sprior'],
                'kindex':aq0para['kindex'],
                'volper':1e-6,
                'onlyafter':aq0para['nrandinit'],
                'check':True,
                'everyn':1
                }
        self.ojfchar = {'dx': len(aq0para['lb']), 'dev': len(aq0para['ev'])}
        self.ojf=f

        self.path = path
        self.fname = fname
        return

C=conf(f,D,n,s,'results','choice.csv')
out = gpbo.search(C)
print out