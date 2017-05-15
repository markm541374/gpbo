import gpbo
import scipy as sp
import copy



#dimension
D=2
#random initialization
ninit=10
#const diagonal noise variance
s=1e-6

def f(x, **ev):
    c=45*sp.exp(-10.*ev['xa'])
    y = -sp.cos(x[0]) - sp.cos(x[1]) + 2.
    b = 0.1*ev['xa'] ** 2
    n = sp.random.normal() * sp.sqrt(s)
    #this is used to check the true value at the imcumbent if available for plotting. Not counted in the optimization
    if 'cheattrue' in ev.keys():
        if ev['cheattrue']:
            b = 0
            n = 0
    print 'f inputs x:{} ev:{} outputs y:{} (b:{} n:{}) c:{}'.format(x, ev, y + n + b, b, n, c)
    return y + b + n, c, dict()


def f_inplane(x,**ev):
    e = copy.deepcopy(ev)
    e['xa']=0
    y, c, aux = f(x, **e)
    return y,c,aux


#pesbs is the predictive entropy search with biased variable cost observations as detailed in https://arxiv.org/abs/1703.04335
C=gpbo.core.config.pesbsdefault(f,D,50,s,'results','pesbsdemo.csv')
C.aqpara['nrandinit']=ninit
C.stoppara = {'tmax': 60 * 60*8}
C.stopfn = gpbo.core.optimize.totaltstopfn
C.aqpara['overhead']='predict'
out = gpbo.search(C)


