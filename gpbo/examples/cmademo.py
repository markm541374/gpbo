import gpbo
import scipy as sp


D=2
n=200
s=0.

def f(x, **ev):
    y = -sp.cos(x[0]-0.05) - sp.cos(x[1]+0.05) + 2
    c = 1.
    n = sp.random.normal() * sp.sqrt(s)
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            n=0
    print 'f inputs x:{} ev:{} outputs y:{} (n:{}) c:{}'.format(x, ev, y + n, n, c)
    return y + n, c, dict()

#import cma
#g = lambda x:-sp.cos(x[0]-0.05) - sp.cos(x[1]+0.05) + 2
#options = {'boundary_handling':'BoundTransform','bounds':[[-1]*D,[1]*D]}
#es = cma.fmin(g,2 * [0], 0.5,options=options)
#print(cma.CMAOptions.defaults())
#raise
class conf():
    """
    fixed s, space is [-1,1]^D

    """
    def __init__(self,f,D,n,s,path,fname):
        self.aqfn = gpbo.core.acquisitions.cmaesaq
        self.aqpara= {
            'ev': {'s': s, 'd': [sp.NaN]},
            'lb': [-1.]*D,
            'ub': [1.]*D
        }

        self.stoppara = {'nmax': n}
        self.stopfn = gpbo.core.optimize.norlocalstopfn

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

C=conf(f,D,n,s,'results','cma.csv')
out = gpbo.search(C)
#print out