# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
from __future__ import print_function
import scipy as sp
from scipy import linalg as spl
import logging
from gpbo.core import GPdc
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

def genmat52ojf(d,lb,ub,ls=0.3):
    from ESutils import gen_dataset
    nt=18
    [X,Y,S,D] = gen_dataset(nt, d, lb, ub, GPdc.MAT52, sp.array([1.5] + [ls] * d))
    G = GPdc.GPcore(X, Y, S, D, GPdc.kernel(GPdc.MAT52, d, sp.array([1.5] + [ls] * d)))
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
    y = spm(spowrap, xmin,  method='l-bfgs-b',bounds=[(-1,1),(-1,1)],options={'ftol':1e-15})
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
    [xmin,ymin,ierror] = DIRECT.solve(dirwrap,lb,ub,user_data=[], algmethod=1, maxf=40000, logfilename='/dev/null')
    #print [xmin, ymin]
    def spowrap(x):
        z = G.infer_m(sp.hstack(sp.array(x) + [0.]), [[sp.NaN]])[0, 0]
        return z
    y = spm(spowrap, xmin, method='l-bfgs-b',bounds=[(-1,1),(-1,1)],options={'ftol':1e-15})
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




def gendecayingpositiveojf(d, lb, ub,sim):

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
        y1=sim*p1(x)+y0
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
