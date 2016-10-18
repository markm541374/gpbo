# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

import scipy as sp
import logging
import GPdc
import DIRECT
from scipy.optimize import minimize as spm
import gpbo
logger = logging.getLogger(__name__)


def trivialojf(x,**ev):
    return sum([(xi-0.1*i)**2 for i,xi in enumerate(x)]),1.,dict()

trivialymin = 0.

def braninojf(x,**ev):
    if 'd' in ev.keys():
        assert(ev['d']==[sp.NaN])
    if 's' in ev.keys():
        noise = sp.random.normal(scale=sp.sqrt(ev['s']))
    else:
        noise=0.
    u = x[0]*7.5 + 2.5
    v = x[1]*7.5 + 2.5
        
    f = (-1.275*(u/sp.pi)**2+5*u/sp.pi+v-6)**2 +(10.-5./(4*sp.pi))*sp.cos(u) + 10.
    
    return f+noise,1.,dict()

braninymin = 0.39788735772973816

def rosenojf(x,**ev):
    if 'd' in ev.keys():
        assert(ev['d']==[sp.NaN])
    if 's' in ev.keys():
        if ev['s']>0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise=0
    else:
        noise=0.
    u = 5.*x[0]
    v = 5.*x[1]
    a=1.
    b=100.
    f = 1e-3*((a-u)**2 + b*(v-u**2)**2)
    return f+noise,1.,dict()

rosenxmin=[0.2,0.2]
rosenymin=0.

def genmat52ojf(d,lb,ub):
    from ESutils import gen_dataset
    nt=18
    [X,Y,S,D] = gen_dataset(nt, d, lb, ub, GPdc.MAT52, sp.array([1.5] + [0.3] * d))
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1.5] + [0.3] * d)))
    def ojf(x,**ev):
        dx=ev['d']
        s=ev['s']
        if ev['s']>0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise=0
        return G.infer_m(sp.array(x),[dx])[0,0]+noise,1.,dict()
    def dirwrap(x,y):
        z = G.infer_m(x,[[sp.NaN]])[0,0]
        #z = obj(x,0.,[sp.NaN])
        return (z,0)
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=20000, logfilename='/dev/null')

    def spowrap(x):
        z = G.infer_m(x,[[sp.NaN]])[0,0]
        #z = obj(x,0.,[sp.NaN])
        return z
    y = spm(spowrap, xmin,  method='nelder-mead',options={'fatol':1e-12})
    xmin = y.x
    ymin = spowrap(y.x)
    logger.info('generated function xmin {} ymin {} globopt:{} locopt:{}'.format(xmin, ymin, ierror, y.status))
    return ojf,xmin,ymin
    
def genbiasedmat52ojf(d,lb,ub,sls):
    #s normalised to 0 exact, 1
    from ESutils import gen_dataset
    nt=20
    [X,Y,S,D] = gen_dataset(nt, d + 1, lb + [0], ub + [1], GPdc.MAT52, sp.array([1.5] + [0.30] * d + [sls]))
    
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d + 1, sp.array([1.5] + [0.30] * d + [sls])))
    def ojf(x,**ev):
        #print "\nojfinput: {} : {}".format(x,ev)
        dx=ev['d']
        s=ev['s']
        if ev['s']>0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise=0
        #print 'noise in ojf {}'.format(noise)
        xa=ev['xa']
        x = sp.array(x)
        xfull = sp.hstack([x,xa])
        return G.infer_m(xfull,[dx])[0,0]+noise,1.,dict()
    
    def dirwrap(x,y):
        z = G.infer_m(sp.hstack(sp.array(x)+[0.]),[[sp.NaN]])[0,0]
        return (z,0)
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=20000, logfilename='/dev/null')
    #print [xmin, ymin]
    def spowrap(x):
        z = G.infer_m(sp.hstack(sp.array(x) + [0.]), [[sp.NaN]])[0, 0]
        return z
    y = spm(spowrap, xmin, method='nelder-mead',options={'fatol':1e-12})
    xmin = y.x
    ymin = spowrap(y.x)
    #print [xmin,ymin]

    if sp.isnan(ymin):
        logger.warning('generator got nan optimizing objective. retrying...')
        return genbiasedmat52ojf(d,lb,ub,sls)
    logger.info('generated function xmin {} ymin {} globopt:{} locopt:{}'.format(xmin, ymin, ierror, y.status))
    return ojf, xmin, ymin

def coswithbias():
    def ojf(x,**ev):
        xa=ev['xa']
        cost = 0.5+xa**2
        y = -sp.cos(x[0]*0.75/sp.pi)*sp.cos(x[1]*0.75/sp.pi)+0.05*xa**2
        return y,cost,dict()
    return ojf,sp.array([0,0]),-1.

def costfnwrap(ojfbase,cfn):
    def ojf(x,**ev):
        y,c0,ojfaux = ojfbase(x,**ev)
        c = cfn(x,**ev)
        return y,c,ojfaux
    return ojf

def cf42(x,**ev):
    return 42.

def cfpower(A,p):
    def cf(x,**ev):
        s=ev['s']
        return A*s**(-p)
    return cf

def cfaexp(A,p):
    def cf(x,**ev):
        s=ev['xa']
        return A*sp.exp(-s*p)
    return cf
class cfnobj():
    def __init__(self,g):
        self.g=g
        return
    def __call__(self,x,**ev):
        xa =  ev['xa']
        return self.g.infer_m(sp.array([[xa]]),[[sp.NaN]])[0,0]

class logcfnobj():
    def __init__(self,g):
        self.g=g
        return
    def __call__(self,x,**ev):
        xa =  ev['xa']
        return sp.exp(self.g.infer_m(sp.array([[xa]]),[[sp.NaN]])[0,0])

def traincfn1d(x,c):
    n = x.size
    g = GPdc.GPcore(x, c, sp.array([1e-1] * n), [[sp.NaN]] * n, GPdc.kernel(GPdc.MAT52, 1, [1., 0.2]))

    if gpbo.core.debugoutput and gpbo.core.debugoptions['cost1d']:
        print 'plotting cost1d...'
        from gpbo.core import debugpath
        import os
        if not os.path.exists(debugpath):
            os.mkdir(debugpath)
        import time
        from matplotlib import pyplot as plt
        f,a=plt.subplots(1)
        low = min(0,min(x))
        high = max(1,max(x))
        xaxis = sp.linspace(low,high,100)
        y,cy = g.infer_diag_post(xaxis,[[sp.NaN]]*100)

        a.plot(xaxis,y[0,:],'b')
        s = 2.*sp.sqrt(cy)
        u=sp.empty(100)
        l=sp.empty(100)
        for i in xrange(100):
            s = sp.sqrt(cy[0,i])
            u[i]=y[0,i]+2.*s
            l[i]=y[0,i]-2.*s
        a.fill_between(xaxis,l,u,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
        for i in xrange(n):
            a.plot(x[i],c[i],'r.')
        f.savefig(os.path.join(debugpath, 'cost1d' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        del(f)
    return cfnobj(g)

def traincfn1dll(x,c):
    #cost modeled in a latent space c=exp(cl)
    n = x.size
    cl=sp.log(c)
    MAP = GPdc.searchMAPhyp(x, cl, sp.array([1e-3] * n), [[sp.NaN]] * n, sp.array([1.,0.,-1.]), sp.array([2.,2.,2.]), GPdc.MAT52CS)
    print 'MAPhyp in costfn {}'.format(MAP)
    g = GPdc.GPcore(x, cl, sp.array([1e-3] * n), [[sp.NaN]] * n, GPdc.kernel(GPdc.MAT52CS,1,MAP))

    if gpbo.core.debugoutput and gpbo.core.debugoptions['cost1d']:
        print 'plotting cost1d...'
        from gpbo.core import debugpath
        import os
        if not os.path.exists(debugpath):
            os.mkdir(debugpath)
        import time
        from matplotlib import pyplot as plt
        f,a=plt.subplots(2)
        low = min(0,min(x))
        high = max(1,max(x))
        xaxis = sp.linspace(low,high,100)
        y,cy = g.infer_diag_post(xaxis,[[sp.NaN]]*100)


        s = 2.*sp.sqrt(cy)
        u=sp.empty(100)
        l=sp.empty(100)
        for i in xrange(100):
            s = sp.sqrt(cy[0,i])
            u[i]=y[0,i]+2.*s
            l[i]=y[0,i]-2.*s
        a[0].plot(xaxis, y[0, :], 'b')
        a[0].fill_between(xaxis,l,u,facecolor='lightblue',edgecolor='lightblue',alpha=0.5)
        for i in xrange(n):
            a[0].plot(x[i],cl[i],'r.')
        a[0].set_ylabel('latent')

        a[1].plot(xaxis, sp.exp(y[0, :]), 'b')
        a[1].fill_between(xaxis, sp.exp(l), sp.exp(u), facecolor='lightblue', edgecolor='lightblue', alpha=0.5)
        for i in xrange(n):
            a[1].plot(x[i], c[i], 'r.')
        a[1].set_ylabel('out')

        f.savefig(os.path.join(debugpath, 'cost1d' + time.strftime('%d_%m_%y_%H:%M:%S') + '.png'))
        del(f)
    return logcfnobj(g)


def gendecayingpositiveojf(d, lb, ub):

    # s normalised to 0 exact, 1
    from ESutils import gen_dataset
    nt = 20
    cl=2.

    [X, Y, S, D] = gen_dataset(nt, d, lb, ub , GPdc.MAT52, sp.array([1] + [0.30] * d))
    G0 = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1] + [0.30] * d)))

    [X, Y, S, D] = gen_dataset(nt, d, lb, ub, GPdc.MAT52, sp.array([1] + [0.30] * d))
    G1 = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1] + [0.30] * d)))

    [X, Y, S, D] = gen_dataset(nt, d, lb, ub, GPdc.MAT52, sp.array([1] + [0.30] * d))
    G2 = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1] + [0.30] * d)))

    def p0(x):
        v= G0.infer_m(x, [[sp.NaN]])[0,0]
        y = v+1 if v>0 else sp.exp(v)
        return y

    def p1(x):
        v= G1.infer_m(x, [[sp.NaN]])[0, 0]
        y = v + 1 if v > 0 else sp.exp(v)
        return y

    def p2(x):
        v= G2.infer_m(x, [[sp.NaN]])[0, 0]
        y = v + 1 if v > 0 else sp.exp(v)
        return y

    def ojf(x, **ev):
        #print "ex: {} {}".format(x,ev['xa'])
        # print "\nojfinput: {} : {}".format(x,ev)
        dx = ev['d']
        s = ev['s']
        if ev['s'] > 0:
            noise = sp.random.normal(scale=sp.sqrt(ev['s']))
        else:
            noise = 0
        # print 'noise in ojf {}'.format(noise)
        xa = ev['xa']
        x = sp.array(x)


        y0=p0(x)
        y1=0.25*p1(x)+y0
        l=p2(x)

        A = (y1-y0)/(sp.exp(l)-1)
        B=y0-A
        y = A*sp.exp(l*xa)+B

        c = sp.exp(-cl * xa)
        return y + noise, c, dict()

    def dirwrap(x, y):
        z = ojf(x,**{'d':[sp.NaN],'s':0,'xa':0})[0]
        return (z, 0)

    [xmin, ymin, ierror] = DIRECT.solve(dirwrap, lb, ub, user_data=[], algmethod=1, maxf=20000, logfilename='/dev/null')

    # print [xmin, ymin]
    def spowrap(x):
        z = ojf(x,**{'d':[sp.NaN],'s':0,'xa':0})[0]
        return z

    y = spm(spowrap, xmin, method='nelder-mead', options={'fatol': 1e-12})
    xmin = y.x
    ymin = spowrap(y.x)
    # print [xmin,ymin]
    logger.info('generated function xmin {} ymin {} yminisnan{} globopt:{} locopt:{}'.format(xmin, ymin, sp.isnan(ymin),ierror, y.status))
    if sp.isnan(ymin):
        logger.warning('generator got nan optimizing objective. retrying...')
        return gendecayingpositiveojf(d, lb, ub)
    return ojf, xmin, ymin
